/**
 * 8 global parameter sliders + seed input.
 * Phase 7, Task P7.3.
 */

export interface GlobalParams {
  seed: number;
  tectonic_activity: number;
  water_abundance: number;
  surface_age: number;
  climate_diversity: number;
  glaciation: number;
  continental_fragmentation: number;
  mountain_prevalence: number;
}

const SLIDER_DEFS = [
  { key: "tectonic_activity",       label: "Tectonic Activity",        min: 0, max: 1, default: 0.35, step: 0.01 },
  { key: "water_abundance",         label: "Water Abundance",           min: 0, max: 1, default: 0.55, step: 0.01 },
  { key: "surface_age",             label: "Surface Age",               min: 0, max: 1, default: 0.50, step: 0.01 },
  { key: "climate_diversity",       label: "Climate Diversity",         min: 0, max: 1, default: 0.70, step: 0.01 },
  { key: "glaciation",              label: "Glaciation",                min: 0, max: 1, default: 0.10, step: 0.01 },
  { key: "continental_fragmentation", label: "Continental Fragmentation", min: 0, max: 1, default: 0.40, step: 0.01 },
  { key: "mountain_prevalence",     label: "Mountain Prevalence",      min: 0, max: 1, default: 0.25, step: 0.01 },
] as const;

export function initSliders(container: HTMLElement): () => GlobalParams {
  const values: Record<string, number> = { seed: 0 };

  // Seed row
  const seedRow = document.createElement("div");
  seedRow.innerHTML = `<label style="font-size:0.8rem">Seed</label><br>
    <input type="number" id="seed-input" value="0" style="width:100%">
    <button id="reroll-btn" style="width:100%;margin-top:4px">Reroll</button>`;
  container.appendChild(seedRow);

  const seedInput = seedRow.querySelector<HTMLInputElement>("#seed-input")!;
  seedRow.querySelector<HTMLButtonElement>("#reroll-btn")!.addEventListener("click", () => {
    const newSeed = Math.floor(Math.random() * 2 ** 32);
    seedInput.value = String(newSeed);
    values["seed"] = newSeed;
  });
  seedInput.addEventListener("change", () => { values["seed"] = Number(seedInput.value); });

  for (const def of SLIDER_DEFS) {
    values[def.key] = def.default;

    const row = document.createElement("div");
    row.style.cssText = "display:flex;flex-direction:column;gap:2px;margin-top:8px";

    const labelEl = document.createElement("label");
    labelEl.style.fontSize = "0.75rem";
    const valueEl = document.createElement("span");
    valueEl.style.cssText = "float:right;font-variant-numeric:tabular-nums";
    valueEl.textContent = def.default.toFixed(2);
    labelEl.append(def.label, " ", valueEl);

    const slider = document.createElement("input");
    slider.type = "range";
    slider.min = String(def.min);
    slider.max = String(def.max);
    slider.step = String(def.step);
    slider.value = String(def.default);
    slider.style.width = "100%";
    slider.addEventListener("input", () => {
      const v = parseFloat(slider.value);
      values[def.key] = v;
      valueEl.textContent = v.toFixed(2);
    });

    row.append(labelEl, slider);
    container.appendChild(row);
  }

  return () => values as unknown as GlobalParams;
}
