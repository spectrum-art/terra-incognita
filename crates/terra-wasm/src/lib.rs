use wasm_bindgen::prelude::*;
use terra_core::generator::{GlobalParams, PlanetGenerator};

#[wasm_bindgen(start)]
pub fn init() {
    // Panic hook for browser console (enabled via optional feature).
    // Add console_error_panic_hook to Cargo.toml when wiring up Phase 7.
}

/// Generate a planet from the given parameters JSON.
/// Returns a flat Float32Array of elevation values + metadata JSON.
/// Phase 7, Task P7.2.
#[wasm_bindgen]
pub async fn generate(params_json: &str) -> Result<JsValue, JsValue> {
    let params: GlobalParams = serde_json::from_str(params_json)
        .map_err(|e| JsValue::from_str(&format!("Invalid params: {e}")))?;

    let generator = PlanetGenerator::new();
    let _result = generator.generate(&params);

    // TODO Phase 7: serialize result to JsValue
    Err(JsValue::from_str("Phase 7: WASM bindings not yet implemented (P7.2)"))
}

/// Compute the realism score for a heightfield.
/// Phase 7, Task P7.2.
#[wasm_bindgen]
pub fn get_score(_heightfield_json: &str) -> Result<JsValue, JsValue> {
    Err(JsValue::from_str("Phase 7: score endpoint not yet implemented (P7.2)"))
}
