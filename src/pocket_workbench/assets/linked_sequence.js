/* Offline sequence and pocket controls for the existing 3Dmol viewer. */
(() => {
  "use strict";

  /** Create a label or control without treating data as markup. */
  function append(parent, tag, text, attributes = {}) {
    const element = document.createElement(tag);
    if (text !== undefined) element.textContent = text;
    for (const [name, value] of Object.entries(attributes)) element.setAttribute(name, value);
    parent.append(element);
    return element;
  }

  /** Compact adjacent positions while preserving gaps in pocket lining. */
  function intervals(positions) {
    const ranges = [];
    for (const position of [...new Set(positions)].sort((a, b) => a - b)) {
      if (ranges.at(-1)?.end === position - 1) ranges.at(-1).end = position;
      else ranges.push({ start: position, end: position });
    }
    return ranges;
  }

  /** Mount the tracks and return the hooks used by the molecular viewer. */
  function mount({ root, model, atoms, onPocket, onResidue, onFeature, onClear }) {
    if (!customElements.get("nightingale-sequence")) return null;
    const axis = model.linked_axis;
    if (!axis) return null;
    const length = axis.sequence.length;
    const byPosition = new Map(), byKey = new Map();
    for (const row of model.residues) {
      const position = row[axis.field], tool = row.tool_residue_number;
      if (!Number.isInteger(position) || position < 1 || position > length || !atoms.has(tool)) {
        continue;
      }
      byPosition.set(position, row);
      byKey.set(row.residue_key, position);
    }
    if (!byPosition.size) return null;

    let start = 1, end = length, picked = null, activeFeature = null;
    let pocketOptions = [], selectedIds = [], rowSignature = "";
    let pocketTrackElements = [];
    const primaryTracks = [], featureButtons = [], actions = new WeakMap();
    const span = append(root, "div", undefined, { class: "sequence-heading" });
    append(span, "h2", "Linked structure and pockets · Nightingale");
    append(span, "span", axis.label, { class: "sequence-badge" });
    const positions = [...byPosition.keys()];
    const first = positions[0], last = positions.at(-1);
    const covered = first === 1 && last === length ? "" :
      ` Coordinates are available for positions ${first}–${last} of ${length}.`;
    append(root, "p", "Choose a lining track, click a residue or enter its position to focus it in 3D. "
      + "Drag the overview window to browse the sequence." + covered, { class: "help" });
    const controls = append(root, "div", undefined, { class: "sequence-controls" });
    const form = append(controls, "form");
    append(form, "label", "Residue position", { for: "linked-residue" });
    const input = append(form, "input", undefined, {
      id: "linked-residue", type: "number", min: 1, max: length, step: 1, required: "",
    });
    append(form, "button", "Find residue", { type: "submit" });
    form.addEventListener("submit", event => {
      event.preventDefault();
      chooseResidue(Number(input.value));
    });
    const full = append(controls, "button", "Full sequence", { type: "button" });
    full.onclick = () => setWindow(1, length);
    const local = append(controls, "button", "Show selected pocket region", { type: "button" });
    local.disabled = true;
    local.onclick = () => {
      const positionList = pocketOptions.filter(p => selectedIds.includes(p.id))
        .flatMap(p => p.positions);
      reveal(positionList);
    };
    const clear = append(controls, "button", "Clear highlights", { type: "button" });
    clear.onclick = () => {
      picked = null;
      input.value = "";
      activeFeature = null;
      message.textContent = "Residue and landmark highlights cleared.";
      updateMarks();
      onClear();
    };
    const message = append(root, "div", "No individual residue selected.", {
      class: "sequence-selection", role: "status", id: "linked-selection",
    });
    const windowLabel = append(root, "div", undefined, { class: "sequence-window" });
    let landmarkColours = { domain: "#6332ac", motif: "#a92370" };

    function makeTrack(parent, name, height, collection) {
      const track = append(parent, `nightingale-${name}`, undefined, {
        length, "display-start": start, "display-end": end, height,
        "highlight-color": "#79c7bd",
      });
      collection.push(track);
      return track;
    }

    function setWindow(a, b) {
      if (!Number.isFinite(a) || !Number.isFinite(b)) return;
      start = Math.max(1, Math.min(Math.max(1, length - 1), Math.floor(a)));
      end = Math.max(Math.min(length, start + 1), Math.min(length, Math.ceil(b)));
      for (const track of [...primaryTracks, ...pocketTrackElements]) {
        if (Number(track.getAttribute("display-start")) !== start) {
          track.setAttribute("display-start", start);
        }
        if (Number(track.getAttribute("display-end")) !== end) {
          track.setAttribute("display-end", end);
        }
      }
      windowLabel.textContent = `Showing ${start}–${end} of ${length} residues`;
    }

    function reveal(positions) {
      if (!positions.length) return;
      let first = length, last = 1;
      for (const position of positions) {
        first = Math.min(first, position);
        last = Math.max(last, position);
      }
      setWindow(first - 12, last + 12);
    }

    function updateMarks() {
      const mark = picked === null ? "" : `${picked}:${picked}`;
      for (const track of [...primaryTracks, ...pocketTrackElements]) {
        if (mark) track.setAttribute("highlight", mark);
        else track.removeAttribute("highlight");
      }
      for (const button of featureButtons) {
        button.element.setAttribute("aria-pressed", button.feature === activeFeature);
      }
    }

    function showResidue(position, zoom = false) {
      if (!Number.isInteger(position) || position < 1 || position > length) return;
      picked = position;
      input.value = String(position);
      const row = byPosition.get(position);
      const author = row?.author_residue_number == null ? "" :
        ` · Original ${row.author_chain || row.tool_chain || ""}:${row.author_residue_number}`
        + (row.insertion_code?.trim() || "");
      const confidence = model.confidence === "none" || row?.plddt == null ? "" :
        ` · pLDDT ${Number(row.plddt).toFixed(1)}`;
      message.textContent = `${axis.sequence[position - 1]}${position}${author}${confidence}`
        + (row ? "" : " · No coordinates in this prepared model");
      if (zoom && (position < start || position > end || end - start > 100)) reveal([position]);
      updateMarks();
    }

    function chooseResidue(position) {
      if (!Number.isInteger(position) || position < 1 || position > length) return;
      showResidue(position, true);
      onResidue(byPosition.get(position)?.tool_residue_number ?? null);
    }

    function pickAxisPosition(event, track) {
      const bounds = track.querySelector("svg")?.getBoundingClientRect();
      if (!bounds?.width) return;
      const fraction = Math.max(0, Math.min(1, (event.clientX - bounds.left) / bounds.width));
      chooseResidue(Math.min(end, start + Math.floor(fraction * (end - start + 1))));
    }

    function selectFeature(feature) {
      activeFeature = activeFeature === feature ? null : feature;
      if (activeFeature) {
        const { start: a, end: b, label } = activeFeature;
        const positions = [...byPosition].filter(([n]) => a <= n && n <= b)
          .map(([, row]) => row.tool_residue_number);
        reveal([a, b]);
        message.textContent = `${label} · ${a}–${b} · ${positions.length}/${b - a + 1} in model`;
        onFeature({ ...activeFeature, toolPositions: positions });
      } else {
        message.textContent = "Landmark highlight cleared.";
        onFeature(null);
      }
      updateMarks();
    }

    function trackFeature({ track, accession, colour, ranges, action }) {
      track.data = [{ accession, color: colour, locations: [{ fragments: ranges }] }];
      actions.set(track, action);
    }

    makeTrack(root, "navigation", 43, primaryTracks);
    const sequence = makeTrack(root, "sequence", 33, primaryTracks);
    sequence.sequence = axis.sequence;
    sequence.addEventListener("click", event => pickAxisPosition(event, sequence));
    if (model.confidence !== "none" && [...byPosition.values()].some(r => r.plddt != null)) {
      const confidence = append(root, "details", undefined, { class: "sequence-confidence" });
      append(confidence, "summary", "Local model confidence (pLDDT, 0–100)");
      const graph = makeTrack(confidence, "linegraph-track", 65, primaryTracks);
      graph.setAttribute("highlight-on-click", "");
      graph.addEventListener("click", event => pickAxisPosition(event, graph));
      graph.data = [{ name: "pLDDT", range: [0, 100], color: "#126c96",
        values: Array.from({ length }, (_, i) => ({ position: i + 1,
          value: byPosition.get(i + 1)?.plddt ?? null })) }];
    }

    if (axis.features.length) {
      append(root, "h3", "Recorded BRCT windows and linker motifs");
      const features = append(root, "div", undefined, { class: "sequence-features" });
      for (const [i, feature] of axis.features.entries()) {
        const row = append(features, "div", undefined, { class: "sequence-track-row" });
        const label = `${feature.label} · ${feature.start}–${feature.end}`;
        const button = append(row, "button", label, { type: "button", "aria-pressed": "false" });
        button.onclick = () => selectFeature(feature);
        const track = makeTrack(row, "track", 18, primaryTracks);
        const updateColour = () => {
          const colour = landmarkColours[feature.kind];
          button.style.setProperty("--feature-colour", colour);
          trackFeature({ track, accession: `landmark-${i}`, colour,
            ranges: [{ start: feature.start, end: feature.end }], action: button.onclick });
        };
        featureButtons.push({ element: button, feature, updateColour });
        updateColour();
      }
    }

    const heading = append(root, "div", undefined, { class: "sequence-heading pocket-track-heading" });
    append(heading, "h3", "Pocket lining tracks");
    const search = append(heading, "input", undefined, {
      type: "search", "aria-label": "Filter pocket tracks by native ID or rank",
      placeholder: "Pocket ID or rank…",
    });
    const rows = append(root, "div", undefined, { class: "sequence-pocket-rows" });
    search.addEventListener("input", drawPocketRows);

    function drawPocketRows() {
      pocketTrackElements = [];
      rows.replaceChildren();
      const term = search.value.trim().toLowerCase();
      const id = term.match(/^(?:pocket\s*)?(\d+)$/)?.[1];
      const shown = pocketOptions.filter(p => !term ||
        (id !== undefined ? p.id === id : p.label.toLowerCase().includes(term)));
      if (!shown.length) {
        append(rows, "p", pocketOptions.length ? "No pocket tracks match this filter." :
          "No pockets are available for the selected method and confidence filters.", { class: "help" });
      }
      if (shown.length > 30) append(rows, "p",
        `Showing the first 30 of ${shown.length} pockets. Use the filter for any pocket ID or rank.`,
        { class: "help" });
      for (const pocket of shown.slice(0, 30)) {
        const row = append(rows, "div", undefined, { class: "sequence-track-row" });
        const button = append(row, "button", pocket.label, {
          type: "button", "data-pocket-id": pocket.id,
          "aria-pressed": String(selectedIds.includes(pocket.id)),
          title: `${pocket.positions.length}/${pocket.total} lining residues mapped to this axis`,
        });
        button.style.setProperty("--feature-colour", pocket.colour);
        const choose = () => {
          reveal(pocket.positions);
          onPocket(pocket.id);
        };
        button.onclick = choose;
        if (pocket.positions.length) {
          trackFeature({ track: makeTrack(row, "track", 22, pocketTrackElements),
            accession: `pocket-${pocket.id}`, colour: pocket.colour,
            ranges: intervals(pocket.positions), action: choose });
        } else {
          append(row, "span", "No mapped lining on this axis", { class: "help" });
        }
      }
      updateMarks();
    }

    root.addEventListener("change", event => {
      const { detail, target } = event;
      if (!detail) return;
      if (detail["display-start"] !== undefined && detail["display-end"] !== undefined) {
        setWindow(Number(detail["display-start"]), Number(detail["display-end"]));
      }
      if ((detail.eventType || detail.eventtype) !== "click") return;
      if (actions.has(target)) actions.get(target)();
      // Axis clicks are handled by the SVG pixel mapping above at all zoom levels.
    });

    root.hidden = false;
    setWindow(1, length);
    return {
      syncPockets(pockets, selected, describe, colour) {
        selectedIds = selected.map(p => String(p.pocket_id));
        pocketOptions = pockets.map(p => {
          const keys = [...new Set(p.residues || [])];
          return { id: String(p.pocket_id), label: describe(p), colour: colour(p),
            total: keys.length, positions: keys.map(key => byKey.get(key)).filter(Number.isInteger) };
        });
        const signature = JSON.stringify(pocketOptions);
        if (signature !== rowSignature) { rowSignature = signature; drawPocketRows(); }
        rows.querySelectorAll("button[data-pocket-id]").forEach(button => {
          button.setAttribute("aria-pressed", selectedIds.includes(button.dataset.pocketId));
        });
        local.disabled = !pocketOptions.some(p => p.positions.length && selectedIds.includes(p.id));
      },
      showModelResidue(key) { showResidue(byKey.get(key), true); },
      focusInitialPocket() { if (!local.disabled) local.click(); },
      setLandmarkColours(colours) {
        landmarkColours = colours;
        for (const button of featureButtons) button.updateColour();
      },
    };
  }

  window.mountPocketWorkbenchSequence = mount;
})();
