/**
 * Realism score panel: 0-100 total + per-metric breakdown.
 * Phase 7, Task P7.5.
 */

export interface MetricScore {
  name: string;
  raw_value: number;
  score_0_1: number;
  passed: boolean;
  subsystem: "noise_synth" | "hydraulic";
}

export interface RealismScoreData {
  total: number;
  metrics: MetricScore[];
}

function totalColor(score: number): string {
  if (score >= 75) return "#55dd66";
  if (score >= 50) return "#ddaa33";
  return "#dd5544";
}

function fmt(v: number, decimals = 3): string {
  return v.toFixed(decimals);
}

export function renderScorePanel(
  panel: HTMLElement,
  data: RealismScoreData,
  generationTimeMs: number,
): void {
  panel.classList.add("visible");

  const totalPct = data.total.toFixed(1);
  const color = totalColor(data.total);

  let html = `
    <div style="font-size:1.2rem;font-weight:700;color:${color};margin-bottom:6px">
      ${totalPct} / 100
    </div>
    <div style="font-size:0.7rem;color:#888;margin-bottom:8px">
      Generated in ${generationTimeMs} ms
    </div>
    <table style="width:100%;border-collapse:collapse;font-size:0.68rem">
      <thead>
        <tr style="color:#aaa;border-bottom:1px solid #333">
          <th style="text-align:left;padding:2px 0">Metric</th>
          <th style="text-align:right;padding:2px 4px">Value</th>
          <th style="text-align:right;padding:2px 0">Score</th>
          <th style="text-align:center;padding:2px 4px">✓</th>
        </tr>
      </thead>
      <tbody>`;

  for (const m of data.metrics) {
    const pct  = (m.score_0_1 * 100).toFixed(0);
    const icon = m.passed ? "✓" : "✗";
    const rowColor = m.passed ? "#e0e0e0" : "#d66";
    const subsysColor = m.subsystem === "noise_synth" ? "#8af" : "#fa8";
    html += `
        <tr style="color:${rowColor};border-bottom:1px solid #1a1a1a">
          <td style="padding:2px 0">
            <span style="display:inline-block;width:6px;height:6px;border-radius:50%;
              background:${subsysColor};margin-right:4px;vertical-align:middle"></span>
            ${m.name}
          </td>
          <td style="text-align:right;padding:2px 4px;font-variant-numeric:tabular-nums">
            ${fmt(m.raw_value)}
          </td>
          <td style="text-align:right;padding:2px 0;font-variant-numeric:tabular-nums">
            ${pct}%
          </td>
          <td style="text-align:center;padding:2px 4px">${icon}</td>
        </tr>`;
  }

  html += `
      </tbody>
    </table>
    <div style="margin-top:6px;font-size:0.62rem;color:#666">
      <span style="color:#8af">&#9679;</span> noise_synth &nbsp;
      <span style="color:#fa8">&#9679;</span> hydraulic
    </div>`;

  panel.innerHTML = html;
}
