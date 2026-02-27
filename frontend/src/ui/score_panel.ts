/**
 * Realism score panel: 0-100 total + per-metric breakdown.
 * Phase 7, Task P7.5.
 */

export interface MetricScore {
  name: string;
  rawValue: number;
  score01: number;
  passed: boolean;
  subsystem: "noise_synth" | "hydraulic";
}

export interface RealismScoreData {
  total: number;
  metrics: MetricScore[];
}

export function renderScorePanel(panel: HTMLElement, data: RealismScoreData): void {
  panel.classList.add("visible");
  panel.innerHTML = `<strong>Realism Score: ${data.total.toFixed(1)}/100</strong><br><hr style="border-color:#444">`;

  for (const m of data.metrics) {
    const colour = m.passed ? "#6d6" : "#d66";
    const pct = (m.score01 * 100).toFixed(0);
    panel.innerHTML += `<div style="color:${colour}">${m.name}: ${pct}%</div>`;
  }
}
