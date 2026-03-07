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
    <div id="pm-header" style="cursor:pointer;display:flex;align-items:center;
      justify-content:space-between;font-weight:700;color:${allColor};font-size:0.9rem;
      user-select:none">
      <span>${allIcon} — ${generationTimeMs} ms</span>
      <span id="pm-chevron" style="font-size:0.7rem;margin-left:6px;
        transition:transform 0.2s;display:inline-block">▼</span>
    </div>
    <div id="pm-collapsible" style="overflow:hidden;transition:max-height 0.25s ease">
    <table style="width:100%;border-collapse:collapse;font-size:0.68rem;margin-top:8px">
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
    </table>
    </div>`; // close #pm-collapsible

  panel.innerHTML = html;

  // Wire collapse toggle
  const header      = panel.querySelector("#pm-header")      as HTMLElement;
  const collapsible = panel.querySelector("#pm-collapsible") as HTMLElement;
  const chevron     = panel.querySelector("#pm-chevron")     as HTMLElement;

  // Set initial max-height (expanded by default)
  collapsible.style.maxHeight = collapsible.scrollHeight + "px";

  header.addEventListener("click", () => {
    const isExpanded = collapsible.style.maxHeight !== "0px";
    collapsible.style.maxHeight = isExpanded ? "0px" : collapsible.scrollHeight + "px";
    chevron.style.transform = isExpanded ? "rotate(-90deg)" : "rotate(0deg)";
  });
}
