//! Hydraulic shaping pipeline: flow routing → stream power → glacial carving
//! → basin delineation.  Phase 6 public API.
pub mod basins;
pub mod flow_routing;
pub mod glacial;
pub mod mass_wasting;
pub mod stream_network;
pub mod stream_power;

use crate::heightfield::HeightField;
use crate::noise::params::{GlacialClass, TerrainClass};
use basins::{delineate_basins, DrainageBasin};
use flow_routing::{compute_d8_flow, FlowField};
use glacial::apply_glacial_carving;
use stream_network::{
    extract_stream_network, StreamNetwork,
    A_MIN_ALPINE, A_MIN_COASTAL, A_MIN_CRATONIC,
    A_MIN_FLUVIAL_ARID, A_MIN_FLUVIAL_HUMID,
};
use stream_power::apply_stream_power;

/// Combined result of one hydraulic shaping pass.
pub struct HydraulicResult {
    pub flow: FlowField,
    pub network: StreamNetwork,
    pub basins: Vec<DrainageBasin>,
}

// ── Per-class parameter tables ────────────────────────────────────────────────

struct HydraulicParams {
    a_min: u32,
    erosion_iters: u32,
    angle_of_repose_deg: f32,
}

fn params_for_class(class: TerrainClass) -> HydraulicParams {
    match class {
        TerrainClass::Alpine => HydraulicParams {
            a_min: A_MIN_ALPINE,
            erosion_iters: 30,
            angle_of_repose_deg: 35.0,
        },
        TerrainClass::FluvialHumid => HydraulicParams {
            a_min: A_MIN_FLUVIAL_HUMID,
            erosion_iters: 50,
            angle_of_repose_deg: 30.0,
        },
        TerrainClass::FluvialArid => HydraulicParams {
            a_min: A_MIN_FLUVIAL_ARID,
            erosion_iters: 20,
            angle_of_repose_deg: 35.0,
        },
        TerrainClass::Cratonic => HydraulicParams {
            a_min: A_MIN_CRATONIC,
            erosion_iters: 10,
            angle_of_repose_deg: 25.0,
        },
        TerrainClass::Coastal => HydraulicParams {
            a_min: A_MIN_COASTAL,
            erosion_iters: 25,
            angle_of_repose_deg: 20.0,
        },
    }
}

// ── Public entry point ────────────────────────────────────────────────────────

/// Apply the full hydraulic shaping pipeline to `hf`.
///
/// Steps:
/// 1. Stream power erosion (iterations and angle of repose are class-specific).
/// 2. Glacial carving (no-op for `GlacialClass::None`).
/// 3. D8 flow routing on the final terrain.
/// 4. Stream network extraction.
/// 5. Drainage basin delineation.
///
/// `erodibility` — per-cell K values in [0, 1]; pass `&[]` for uniform K=0.5.
pub fn apply_hydraulic_shaping(
    hf: &mut HeightField,
    terrain_class: TerrainClass,
    erodibility: &[f32],
    glacial_class: GlacialClass,
) -> HydraulicResult {
    let p = params_for_class(terrain_class);

    // Step 1 — stream power erosion.  Returns the final flow field after the
    // last erosion iteration.
    let flow_after_erosion = apply_stream_power(hf, erodibility, p.erosion_iters, p.angle_of_repose_deg);

    // Step 2 — glacial carving (borrows pre-erosion flow field only for the
    // glacial mask; recomputes internally after carving).
    apply_glacial_carving(hf, &flow_after_erosion, glacial_class);

    // Step 3 — final flow routing on the shaped terrain.
    let flow = compute_d8_flow(hf);

    // Step 4 — stream network.
    let network = extract_stream_network(&flow, p.a_min);

    // Step 5 — basin delineation.
    let basins = delineate_basins(&flow, hf);

    HydraulicResult { flow, network, basins }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::heightfield::HeightField;

    fn make_ramp(rows: usize, cols: usize) -> HeightField {
        let deg = cols as f64 * 0.0009;
        let mut hf = HeightField::new(cols, rows, 0.0, deg, 0.0, deg, 0.0);
        for r in 0..rows {
            for c in 0..cols {
                hf.set(r, c, (cols - c) as f32 * 20.0);
            }
        }
        hf
    }

    /// V-valley: walls slope toward a central channel that drains southward.
    /// Centre column accumulation ≈ rows × (cols/2) — well above any A_min.
    fn make_valley(rows: usize, cols: usize) -> HeightField {
        let center = cols / 2;
        let deg = cols as f64 * 0.0009;
        let mut hf = HeightField::new(cols, rows, 0.0, deg, 0.0, deg, 0.0);
        for r in 0..rows {
            for c in 0..cols {
                let dist = (c as isize - center as isize).unsigned_abs() as f32;
                hf.set(r, c, dist * 100.0 + (rows - 1 - r) as f32 * 50.0 + 1000.0);
            }
        }
        hf
    }

    #[test]
    fn basin_areas_cover_all_cells() {
        let mut hf = make_ramp(16, 32);
        let result = apply_hydraulic_shaping(
            &mut hf,
            TerrainClass::FluvialHumid,
            &[],
            GlacialClass::None,
        );
        let total: u32 = result.basins.iter().map(|b| b.area_cells).sum();
        assert_eq!(
            total,
            (16 * 32) as u32,
            "basin areas must sum to total cells, got {total}"
        );
    }

    #[test]
    fn stream_network_non_empty_after_shaping() {
        // V-valley: cells from both sides converge to the centre column, giving
        // accumulation ≈ rows × (cols/2) >> A_min for all terrain classes.
        let mut hf = make_valley(32, 32);
        let result = apply_hydraulic_shaping(
            &mut hf,
            TerrainClass::FluvialHumid,
            &[],
            GlacialClass::None,
        );
        assert!(
            result.network.max_order >= 1,
            "stream network must have at least Strahler order 1"
        );
        let stream_count = result.network.stream_cells.iter().filter(|&&s| s).count();
        assert!(stream_count > 0, "no stream cells found after shaping");
    }

    #[test]
    fn all_terrain_classes_complete_without_panic() {
        use TerrainClass::*;
        for tc in [Alpine, FluvialHumid, FluvialArid, Cratonic, Coastal] {
            let mut hf = make_ramp(8, 16);
            let result = apply_hydraulic_shaping(&mut hf, tc, &[], GlacialClass::None);
            let total: u32 = result.basins.iter().map(|b| b.area_cells).sum();
            assert_eq!(total, (8 * 16) as u32, "class {tc:?}: basin sum mismatch");
        }
    }
}
