//! WASM bindings for the terra-core generation pipeline.
//! Phase 7, Task P7.2.

use serde::{Deserialize, Serialize};
use terra_core::generator::{
    derive_debug_params, generate_at_location, GlobalParams, PlanetGenerator,
};
use terra_core::metrics::score::RealismScore;
use terra_core::noise::params::GlacialClass;
use terra_core::planet::{generate_planet_overview, OVERVIEW_HEIGHT, OVERVIEW_WIDTH};
use terra_core::plates::regime_field::TectonicRegime;
use wasm_bindgen::prelude::*;

#[wasm_bindgen(start)]
pub fn init() {
    // Panic hook wired when console_error_panic_hook is added (Phase 8).
}

// ── Serialisable mirror types ─────────────────────────────────────────────────

#[derive(Serialize, Deserialize)]
struct MetricScoreJs {
    name: String,
    raw_value: f32,
    score_0_1: f32,
    passed: bool,
    subsystem: String,
}

#[derive(Serialize)]
struct RealismScoreJs {
    total: f32,
    metrics: Vec<MetricScoreJs>,
}

#[derive(Serialize)]
struct PlanetResultJs {
    heights: Vec<f32>,
    regimes: Vec<u8>,
    map_field: Vec<f32>,
    width: u32,
    height: u32,
    score: RealismScoreJs,
    generation_time_ms: u64,
}

fn regime_to_u8(r: TectonicRegime) -> u8 {
    match r {
        TectonicRegime::PassiveMargin => 0,
        TectonicRegime::CratonicShield => 1,
        TectonicRegime::ActiveCompressional => 2,
        TectonicRegime::ActiveExtensional => 3,
        TectonicRegime::VolcanicHotspot => 4,
    }
}

fn score_to_js(s: RealismScore) -> RealismScoreJs {
    RealismScoreJs {
        total: s.total,
        metrics: s
            .metrics
            .into_iter()
            .map(|m| MetricScoreJs {
                name: m.name.to_owned(),
                raw_value: m.raw_value,
                score_0_1: m.score_0_1,
                passed: m.passed,
                subsystem: m.subsystem.to_owned(),
            })
            .collect(),
    }
}

// ── Public WASM API ───────────────────────────────────────────────────────────

#[wasm_bindgen]
pub fn generate(params_js: JsValue) -> Result<JsValue, JsValue> {
    let params: GlobalParams = serde_wasm_bindgen::from_value(params_js)
        .map_err(|e| JsValue::from_str(&format!("Invalid params: {e}")))?;

    // Use JS Date.now() for timing — std::time::Instant panics on wasm32.
    let t0 = js_sys::Date::now();
    let result = PlanetGenerator::new().generate(&params);
    let generation_time_ms = (js_sys::Date::now() - t0) as u64;

    let js_result = PlanetResultJs {
        heights: result.heightfield.data,
        regimes: result.regime_field.into_iter().map(regime_to_u8).collect(),
        map_field: result.map_field,
        width: terra_core::generator::GRID_WIDTH as u32,
        height: terra_core::generator::GRID_HEIGHT as u32,
        score: score_to_js(result.score),
        generation_time_ms,
    };

    serde_wasm_bindgen::to_value(&js_result)
        .map_err(|e| JsValue::from_str(&format!("Serialisation error: {e}")))
}

// ── Planet overview types ─────────────────────────────────────────────────────

#[derive(Serialize)]
struct PlanetMetricJs {
    name: String,
    raw_value: f32,
    threshold: f32,
    pass: bool,
    description: String,
}

#[derive(Serialize)]
struct PlanetMetricsJs {
    metrics: Vec<PlanetMetricJs>,
    all_pass: bool,
}

#[derive(Serialize)]
struct PlanetOverviewJs {
    elevations: Vec<f32>,
    /// true = ocean pixel.
    ocean_mask: Vec<bool>,
    sea_level_km: f32,
    regimes: Vec<u8>,
    map_field: Vec<f32>,
    erodibility_field: Vec<f32>,
    /// 0 = None, 1 = Former, 2 = Active.
    glaciation: Vec<u8>,
    planet_metrics: PlanetMetricsJs,
    width: u32,
    height: u32,
    generation_time_ms: u64,
}

fn glacial_to_u8(g: GlacialClass) -> u8 {
    match g {
        GlacialClass::None => 0,
        GlacialClass::Former => 1,
        GlacialClass::Active => 2,
    }
}

// ── generate_overview binding ─────────────────────────────────────────────────

/// Generate a full 1024×512 planet overview from `GlobalParams`.
///
/// Separate from `generate()` — the existing tile pipeline is unchanged.
#[wasm_bindgen]
pub fn generate_overview(params_js: JsValue) -> Result<JsValue, JsValue> {
    let params: GlobalParams = serde_wasm_bindgen::from_value(params_js)
        .map_err(|e| JsValue::from_str(&format!("Invalid params: {e}")))?;

    let t0 = js_sys::Date::now();
    let overview = generate_planet_overview(&params);
    let generation_time_ms = (js_sys::Date::now() - t0) as u64;

    let js = PlanetOverviewJs {
        elevations: overview.elevations,
        ocean_mask: overview.ocean_mask,
        sea_level_km: overview.sea_level_km,
        regimes: overview.regimes.into_iter().map(regime_to_u8).collect(),
        map_field: overview.map_field,
        erodibility_field: overview.erodibility_field,
        glaciation: overview.glaciation.into_iter().map(glacial_to_u8).collect(),
        planet_metrics: PlanetMetricsJs {
            all_pass: overview.planet_metrics.all_pass,
            metrics: overview
                .planet_metrics
                .metrics
                .iter()
                .map(|m| PlanetMetricJs {
                    name: m.name.to_owned(),
                    raw_value: m.raw_value,
                    threshold: m.threshold,
                    pass: m.pass,
                    description: m.description.to_owned(),
                })
                .collect(),
        },
        width: OVERVIEW_WIDTH as u32,
        height: OVERVIEW_HEIGHT as u32,
        generation_time_ms,
    };

    serde_wasm_bindgen::to_value(&js)
        .map_err(|e| JsValue::from_str(&format!("Serialisation error: {e}")))
}

// ── generate_at_location binding ──────────────────────────────────────────────

#[derive(Serialize)]
struct SampledFieldsJs {
    terrain_class: String,
    local_regime: String,
    local_map_mm: f32,
    local_erodibility: f32,
    local_grain_angle: f32,
    local_grain_intensity: f32,
    local_glaciation: String,
}

#[derive(Serialize)]
struct LocationTileResultJs {
    heights: Vec<f32>,
    regimes: Vec<u8>,
    map_field: Vec<f32>,
    width: u32,
    height: u32,
    score: RealismScoreJs,
    generation_time_ms: u64,
    lat: f32,
    lon: f32,
    sampled_fields: SampledFieldsJs,
}

fn terrain_class_to_str(tc: terra_core::noise::params::TerrainClass) -> String {
    format!("{tc:?}")
}

fn regime_to_str(r: TectonicRegime) -> String {
    match r {
        TectonicRegime::PassiveMargin => "PassiveMargin",
        TectonicRegime::CratonicShield => "CratonicShield",
        TectonicRegime::ActiveCompressional => "ActiveCompressional",
        TectonicRegime::ActiveExtensional => "ActiveExtensional",
        TectonicRegime::VolcanicHotspot => "VolcanicHotspot",
    }
    .to_owned()
}

fn glacial_to_str(g: GlacialClass) -> String {
    match g {
        GlacialClass::None => "None",
        GlacialClass::Former => "Former",
        GlacialClass::Active => "Active",
    }
    .to_owned()
}

/// Generate a tile characterised by the planet fields at the given lat/lon.
///
/// Runs the full planet simulation at 1024×512, samples spatial fields at the
/// clicked cell, then runs the tile pipeline at 512×256.
#[wasm_bindgen]
pub fn generate_at_location_wasm(
    params_js: JsValue,
    lat: f32,
    lon: f32,
) -> Result<JsValue, JsValue> {
    let params: GlobalParams = serde_wasm_bindgen::from_value(params_js)
        .map_err(|e| JsValue::from_str(&format!("Invalid params: {e}")))?;

    let t0 = js_sys::Date::now();
    let r = generate_at_location(&params, lat, lon);
    let generation_time_ms = (js_sys::Date::now() - t0) as u64;

    let js = LocationTileResultJs {
        heights: r.heightfield.data,
        regimes: r.regime_field.into_iter().map(regime_to_u8).collect(),
        map_field: r.map_field,
        width: terra_core::generator::GRID_WIDTH as u32,
        height: terra_core::generator::GRID_HEIGHT as u32,
        score: score_to_js(r.score),
        generation_time_ms,
        lat: r.lat,
        lon: r.lon,
        sampled_fields: SampledFieldsJs {
            terrain_class: terrain_class_to_str(r.terrain_class),
            local_regime: regime_to_str(r.local_regime),
            local_map_mm: r.local_map_mm,
            local_erodibility: r.local_erodibility,
            local_grain_angle: r.local_grain_angle,
            local_grain_intensity: r.local_grain_intensity,
            local_glaciation: glacial_to_str(r.local_glaciation),
        },
    };

    serde_wasm_bindgen::to_value(&js)
        .map_err(|e| JsValue::from_str(&format!("Serialisation error: {e}")))
}

/// Resolve GlobalParams → internal DebugParams without running the full pipeline.
///
/// Use this to verify slider wiring: each slider should change at least one
/// field when moved from 0 to 1. Prints-to-console from JS:
/// ```js
/// console.log(debug_params(params));
/// ```
#[wasm_bindgen]
pub fn debug_params(params_js: JsValue) -> Result<JsValue, JsValue> {
    let params: GlobalParams = serde_wasm_bindgen::from_value(params_js)
        .map_err(|e| JsValue::from_str(&format!("Invalid params: {e}")))?;

    let dbg = derive_debug_params(&params);
    serde_wasm_bindgen::to_value(&dbg)
        .map_err(|e| JsValue::from_str(&format!("Serialisation error: {e}")))
}
