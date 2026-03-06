/**
 * Planet-scale metrics panel: 6-metric spatial test battery.
 * Phase A, PA.4 display.
 */

export interface PlanetMetric {
  name:        string;
  raw_value:   number;
  threshold:   number;
  pass:        boolean;
  description: string;
}

export interface PlanetMetricsData {
  metrics:  PlanetMetric[];
  all_pass: boolean;
}

function fmt(v: number): string {
  return v.toFixed(3);
}

export function renderPlanetMetricsPanel(
  panel: HTMLElement,
  data:  PlanetMetricsData,
  generationTimeMs: number,
): void {
  panel.classList.add("visible");

  const allColor = data.all_pass ? "#55dd66" : "#dd9933";
  const allIcon  = data.all_pass ? "✓ All pass" : "⚠ Some metrics failed";

  let html = `
    <div style="font-weight:700;color:${allColor};margin-bottom:6px;font-size:0.9rem">
      ${allIcon}
    </div>
    <div style="font-size:0.7rem;color:#888;margin-bottom:8px">
      Planet overview — ${generationTimeMs} ms
    </div>
    <table style="width:100%;border-collapse:collapse;font-size:0.68rem">
      <thead>
        <tr style="color:#aaa;border-bottom:1px solid #333">
          <th style="text-align:left;padding:2px 0">Metric</th>
          <th style="text-align:right;padding:2px 4px">Value</th>
          <th style="text-align:center;padding:2px 4px">✓</th>
        </tr>
      </thead>
      <tbody>`;

  for (const m of data.metrics) {
    const icon     = m.pass ? "✓" : "✗";
    const rowColor = m.pass ? "#e0e0e0" : "#d66";
    html += `
        <tr style="color:${rowColor};border-bottom:1px solid #1a1a1a"
            title="${m.description}">
          <td style="padding:2px 0">${m.name}</td>
          <td style="text-align:right;padding:2px 4px;font-variant-numeric:tabular-nums">
            ${fmt(m.raw_value)}
          </td>
          <td style="text-align:center;padding:2px 4px">${icon}</td>
        </tr>`;
  }

  html += `
      </tbody>
    </table>`;

  panel.innerHTML = html;
}
