//! WASM bindings for the terra-core generation pipeline.
//! Phase 7, Task P7.2.

use wasm_bindgen::prelude::*;
use serde::{Deserialize, Serialize};
use terra_core::generator::{GlobalParams, PlanetGenerator};
use terra_core::metrics::score::{compute_realism_score, RealismScore};
use terra_core::noise::params::TerrainClass;
use terra_core::plates::regime_field::TectonicRegime;

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
        TectonicRegime::PassiveMargin       => 0,
        TectonicRegime::CratonicShield      => 1,
        TectonicRegime::ActiveCompressional => 2,
        TectonicRegime::ActiveExtensional   => 3,
        TectonicRegime::VolcanicHotspot     => 4,
    }
}

fn score_to_js(s: RealismScore) -> RealismScoreJs {
    RealismScoreJs {
        total: s.total,
        metrics: s.metrics.into_iter().map(|m| MetricScoreJs {
            name: m.name.to_owned(),
            raw_value: m.raw_value,
            score_0_1: m.score_0_1,
            passed: m.passed,
            subsystem: m.subsystem.to_owned(),
        }).collect(),
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

#[wasm_bindgen]
pub fn get_score(heightfield_js: JsValue) -> Result<JsValue, JsValue> {
    #[derive(Deserialize)]
    struct Input {
        heights: Vec<f32>,
        width: u32,
        height: u32,
        terrain_class: String,
    }

    let input: Input = serde_wasm_bindgen::from_value(heightfield_js)
        .map_err(|e| JsValue::from_str(&format!("Invalid input: {e}")))?;

    let tc = match input.terrain_class.as_str() {
        "Alpine"       => TerrainClass::Alpine,
        "FluvialArid"  => TerrainClass::FluvialArid,
        "Cratonic"     => TerrainClass::Cratonic,
        "Coastal"      => TerrainClass::Coastal,
        _              => TerrainClass::FluvialHumid,
    };

    let hf = terra_core::heightfield::HeightField {
        data: input.heights,
        width: input.width as usize,
        height: input.height as usize,
        min_lon: -180.0,
        max_lon:  180.0,
        min_lat:  -90.0,
        max_lat:   90.0,
    };

    let score = compute_realism_score(&hf, tc);
    serde_wasm_bindgen::to_value(&score_to_js(score))
        .map_err(|e| JsValue::from_str(&format!("Serialisation error: {e}")))
}
