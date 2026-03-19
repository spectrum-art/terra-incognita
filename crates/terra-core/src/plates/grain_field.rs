//! Structural grain field derived from continuous boundary character.
//!
//! Note: arc/ridge geometry grain derivation was removed in the plate system rebuild.
//! See git history before commit `eb343e4` for the previous implementation.

use crate::plates::age_field::{cell_to_vec3, distance_to_seeds_km};
use crate::plates::plate_dynamics::PlateDynamics;
use crate::plates::regime_field::{RegimeCharacterField, RegimeField, TectonicRegime};
use crate::sphere::Vec3;

/// Structural grain vector field: angle (radians) and intensity [0, 1] per cell.
#[derive(Clone, Debug, PartialEq)]
pub struct GrainField {
    pub angles: Vec<f32>,
    pub intensities: Vec<f32>,
    pub width: usize,
    pub height: usize,
}

impl GrainField {
    pub fn zero(width: usize, height: usize) -> Self {
        Self {
            angles: vec![0.0; width * height],
            intensities: vec![0.0; width * height],
            width,
            height,
        }
    }
}

pub fn derive_grain_field(
    regime_character: &RegimeCharacterField,
    regime_field: &RegimeField,
    dynamics: &PlateDynamics,
    hotspots: &[Vec3],
    width: usize,
    height: usize,
) -> GrainField {
    let mut field = GrainField::zero(width, height);
    let boundary_seeds: Vec<usize> = dynamics
        .is_boundary
        .iter()
        .enumerate()
        .filter_map(|(idx, &is_boundary)| is_boundary.then_some(idx))
        .collect();
    let nearest_boundary = distance_to_seeds_km(width, height, &boundary_seeds);

    for row in 0..height {
        for col in 0..width {
            let idx = row * width + col;
            if regime_field.data[idx] == TectonicRegime::CratonicShield {
                continue;
            }

            let dominant = dominant_boundary_mode(regime_character, idx);
            let point = cell_to_vec3(row, col, width, height);
            let (angle, intensity) = match dominant {
                DominantMode::Convergent | DominantMode::Divergent | DominantMode::Transform => {
                    let source = nearest_boundary.nearest_source[idx];
                    if source == usize::MAX {
                        (0.0, 0.0)
                    } else {
                        let tangent = (
                            dynamics.boundary_field[source].tangent_east,
                            dynamics.boundary_field[source].tangent_north,
                        );
                        let tangent_angle = tangent.1.atan2(tangent.0);
                        let angle = match dominant {
                            DominantMode::Convergent => tangent_angle + std::f32::consts::FRAC_PI_2,
                            DominantMode::Divergent | DominantMode::Transform => tangent_angle,
                            DominantMode::Hotspot => unreachable!(),
                        };
                        let intensity = match dominant {
                            DominantMode::Convergent => regime_character.convergent_influence[idx],
                            DominantMode::Divergent => regime_character.divergent_influence[idx],
                            DominantMode::Transform => regime_character.transform_influence[idx],
                            DominantMode::Hotspot => 0.0,
                        };
                        (angle, intensity)
                    }
                }
                DominantMode::Hotspot => {
                    let (angle, intensity) = nearest_hotspot_angle(point, hotspots);
                    (angle, intensity * regime_character.hotspot_influence[idx])
                }
            };
            field.angles[idx] = angle;
            field.intensities[idx] = intensity.clamp(0.0, 1.0);
        }
    }

    field
}

#[derive(Clone, Copy)]
enum DominantMode {
    Convergent,
    Divergent,
    Transform,
    Hotspot,
}

fn dominant_boundary_mode(character: &RegimeCharacterField, idx: usize) -> DominantMode {
    let convergent = character.convergent_influence[idx];
    let divergent = character.divergent_influence[idx];
    let transform = character.transform_influence[idx];
    let hotspot = character.hotspot_influence[idx];
    if convergent >= divergent && convergent >= transform && convergent >= hotspot {
        DominantMode::Convergent
    } else if divergent >= transform && divergent >= hotspot {
        DominantMode::Divergent
    } else if transform >= hotspot {
        DominantMode::Transform
    } else {
        DominantMode::Hotspot
    }
}

fn nearest_hotspot_angle(point: Vec3, hotspots: &[Vec3]) -> (f32, f32) {
    let Some((&hotspot, distance)) = hotspots
        .iter()
        .map(|hotspot| (hotspot, point.dot(*hotspot).clamp(-1.0, 1.0).acos() as f32))
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
    else {
        return (0.0, 0.0);
    };

    let (lat_deg, lon_deg) = point.to_latlon();
    let lat = lat_deg.to_radians();
    let lon = lon_deg.to_radians();
    let east = Vec3::new(-lon.sin(), lon.cos(), 0.0);
    let north = Vec3::new(-lat.sin() * lon.cos(), -lat.sin() * lon.sin(), lat.cos());
    let tangent = Vec3::new(
        hotspot.x - point.x * point.dot(hotspot),
        hotspot.y - point.y * point.dot(hotspot),
        hotspot.z - point.z * point.dot(hotspot),
    );
    let angle = tangent.dot(north).atan2(tangent.dot(east)) as f32;
    let intensity = (1.0 - distance / (300.0_f32 / 6371.0_f32)).clamp(0.0, 1.0);
    (angle, intensity)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plates::plate_dynamics::{BoundaryCharacter, PlateDynamics};

    fn sample_dynamics() -> PlateDynamics {
        let mut boundary_field = vec![BoundaryCharacter::default(); 4];
        boundary_field[1].tangent_east = 1.0;
        boundary_field[1].tangent_north = 0.0;
        boundary_field[1].convergent_rate = 4.0;
        boundary_field[2].tangent_east = 1.0;
        boundary_field[2].tangent_north = 0.0;
        boundary_field[2].convergent_rate = -4.0;
        PlateDynamics {
            plate_velocities: vec![(0.0, 0.0); 2],
            boundary_field,
            is_boundary: vec![false, true, true, false],
        }
    }

    #[test]
    fn cratons_have_zero_intensity() {
        let character = RegimeCharacterField {
            convergent_influence: vec![0.0; 4],
            divergent_influence: vec![0.0; 4],
            transform_influence: vec![0.0; 4],
            hotspot_influence: vec![0.0; 4],
            cratonic_stability: vec![1.0; 4],
            width: 4,
            height: 1,
        };
        let regime = RegimeField {
            data: vec![TectonicRegime::CratonicShield; 4],
            width: 4,
            height: 1,
        };
        let grain = derive_grain_field(&character, &regime, &sample_dynamics(), &[], 4, 1);
        assert!(grain.intensities.iter().all(|&value| value == 0.0));
    }

    #[test]
    fn convergent_grain_is_perpendicular_to_tangent() {
        let character = RegimeCharacterField {
            convergent_influence: vec![0.0, 1.0, 0.0, 0.0],
            divergent_influence: vec![0.0; 4],
            transform_influence: vec![0.0; 4],
            hotspot_influence: vec![0.0; 4],
            cratonic_stability: vec![0.0; 4],
            width: 4,
            height: 1,
        };
        let regime = RegimeField {
            data: vec![TectonicRegime::PassiveMargin; 4],
            width: 4,
            height: 1,
        };
        let grain = derive_grain_field(&character, &regime, &sample_dynamics(), &[], 4, 1);
        let direction = (grain.angles[1].cos(), grain.angles[1].sin());
        let tangent = (1.0_f32, 0.0_f32);
        let dot = direction.0 * tangent.0 + direction.1 * tangent.1;
        assert!(dot.abs() < 0.2);
    }

    #[test]
    fn divergent_grain_is_parallel_to_tangent() {
        let character = RegimeCharacterField {
            convergent_influence: vec![0.0; 4],
            divergent_influence: vec![0.0, 0.0, 1.0, 0.0],
            transform_influence: vec![0.0; 4],
            hotspot_influence: vec![0.0; 4],
            cratonic_stability: vec![0.0; 4],
            width: 4,
            height: 1,
        };
        let regime = RegimeField {
            data: vec![TectonicRegime::PassiveMargin; 4],
            width: 4,
            height: 1,
        };
        let grain = derive_grain_field(&character, &regime, &sample_dynamics(), &[], 4, 1);
        let direction = (grain.angles[2].cos(), grain.angles[2].sin());
        let tangent = (1.0_f32, 0.0_f32);
        let dot = direction.0 * tangent.0 + direction.1 * tangent.1;
        assert!(dot.abs() > 0.8);
    }
}
