/* Literal protein sequence search and independent 3D highlight controls. */
((global) => {
  "use strict";

  /** Find every exact, one-based occurrence, including overlapping occurrences. */
  function motifHits(sequence, input) {
    const clean = input.replace(/\s/g, "");
    if (!/^[A-Za-z]{1,100}$/.test(clean)) {
      throw new Error("Enter 1–100 amino acid letters. Search is exact; X is a literal letter.");
    }
    const motif = clean.toUpperCase();
    const hits = [];
    for (let i = sequence.indexOf(motif); i !== -1; i = sequence.indexOf(motif, i + 1)) {
      hits.push({ start: i + 1, end: i + motif.length });
    }
    return { motif, hits };
  }

  /** Map inclusive sequence spans to the verified residues present in the viewer. */
  function mappedResidues(rows, field, spans) {
    const selected = new Set();
    const sorted = [...spans].sort((a, b) => a.start - b.start);
    let i = 0, through = 0;
    for (const row of rows) {
      const position = row[field];
      if (!Number.isInteger(position) || !Number.isInteger(row.tool_residue_number)) continue;
      while (i < sorted.length && sorted[i].start <= position) {
        through = Math.max(through, sorted[i].end);
        i += 1;
      }
      if (position <= through) selected.add(row.tool_residue_number);
    }
    return [...selected];
  }

  /** Make safe DOM controls without parsing imported sequence/model text as HTML. */
  function append(parent, tag, text, attrs = {}) {
    const element = document.createElement(tag);
    if (text) element.textContent = text;
    for (const [key, value] of Object.entries(attrs)) element.setAttribute(key, value);
    parent.append(element);
    return element;
  }

  /** Add one search, hit selector and colour picker to the shared structure viewer. */
  function mount({ root, models, states, onMotifs, onLandmark, onLandmarkColour }) {
    if (!models.some(model => model.linked_axis)) return null;
    root.hidden = false;
    const header = append(root, "h2", "Amino acid search");
    header.id = "motif-title";
    append(root, "p", "Search exact amino acid letters and highlight one occurrence or all of them. "
      + "Pocket and domain highlights stay selected.", { class: "help" });
    const form = append(root, "form", "", { class: "motif-form" });
    const target = append(form, "select", "", { "aria-label": "Model for amino acid search" });
    target.hidden = models.length === 1;
    const available = models.flatMap((model, i) => model.linked_axis ? [i] : []);
    if (available.length > 1) append(target, "option", "Both models", { value: "both" });
    models.forEach((model, i) => {
      const option = append(target, "option", model.short_label || model.label, { value: i });
      option.disabled = !model.linked_axis;
    });
    target.value = available.length > 1 ? "both" : String(available[0]);
    const input = append(form, "input", "", { "aria-label": "Amino acid sequence", type: "text",
      autocomplete: "off", spellcheck: "false", maxlength: 200, placeholder: "e.g. ESQDR" });
    append(form, "button", "Find amino acids", { type: "submit" });
    const colourLabel = append(form, "label", "Motif search colour");
    const colour = append(colourLabel, "input", "", { type: "color", value: "#bd315b" });
    const select = append(form, "select", "", { "aria-label": "Amino acid occurrence" });
    select.hidden = true;
    const clear = append(form, "button", "Clear motif search", { type: "button" });
    clear.disabled = true;
    const status = append(root, "p", "", { class: "motif-result", role: "status", id: "motif-result" });
    let search = null;
    const scopes = () => target.value === "both" ? available : [Number(target.value)];

    function selection({ reframe }) {
      if (!search?.hits.length) { onMotifs([], false); return; }
      const all = select.value === "all";
      const chosen = all ? search.hits : [search.hits[Number(select.value)]];
      const selections = scopes().flatMap(i => {
        const spans = chosen.filter(hit => hit.index === i);
        if (!spans.length) return [];
        const state = states[i], axis = state.model.linked_axis;
        const toolPositions = mappedResidues(state.model.residues, axis.field, spans)
          .filter(position => state.atoms.has(position));
        const label = !all || spans.length === 1 ? `${search.motif} · ${spans[0].start}–${spans[0].end}` :
          `${search.motif} · ${spans.length} occurrences`;
        return [{ index: i, spans, label, colour: colour.value, toolPositions }];
      });
      const hitLabel = search.hits.length === 1 ? "1 exact occurrence" : `${search.hits.length} exact occurrences`;
      const details = selections.map(({ index, label, toolPositions }) => {
        const model = models[index];
        const prefix = models.length > 1 ? `${model.short_label} · ` : "";
        const mapped = toolPositions.length ? `${toolPositions.length} residues with coordinates` :
          "no coordinates for this selection";
        return `${prefix}${label} · ${mapped} · ${model.linked_axis.label}`;
      });
      status.textContent = `${hitLabel} · ${details.join(" | ")}`;
      if (search.hits.length > 1000) status.textContent += " · first 1,000 individually selectable";
      onMotifs(selections, reframe);
    }

    function runSearch() {
      const indices = scopes();
      select.replaceChildren();
      select.hidden = true;
      clear.disabled = !input.value;
      search = null;
      if (!input.value.trim()) { status.textContent = "Enter an amino acid string."; onMotifs([], false); return; }
      try {
        const searches = indices.map(index => ({ index, ...motifHits(models[index].linked_axis.sequence, input.value) }));
        search = { motif: searches[0].motif,
          hits: searches.flatMap(({ index, hits }) => hits.map(hit => ({ ...hit, index }))) };
      } catch (error) { status.textContent = error.message; onMotifs([], false); return; }
      if (!search.hits.length) {
        status.textContent = `No exact occurrence of ${search.motif} in the selected sequence${indices.length > 1 ? "s" : ""}.`;
        onMotifs([], false);
        return;
      }
      append(select, "option", "All occurrences", { value: "all" });
      search.hits.slice(0, 1000).forEach(({ index, start, end }, i) => {
        const prefix = indices.length > 1 ? `${models[index].short_label} · ` : "";
        append(select, "option", `${prefix}${i + 1} · ${start}–${end}`, { value: i });
      });
      select.hidden = false;
      selection({ reframe: true });
    }

    form.addEventListener("submit", event => { event.preventDefault(); runSearch(); });
    target.addEventListener("change", () => { if (input.value) runSearch(); });
    select.addEventListener("change", () => selection({ reframe: true }));
    colour.addEventListener("change", () => selection({ reframe: false }));
    clear.addEventListener("click", () => {
      input.value = ""; search = null; select.replaceChildren(); select.hidden = true;
      status.textContent = "Motif search cleared."; clear.disabled = true; onMotifs([], false);
    });

    // The single-model Nightingale view already supplies clickable domain tracks.
    if (models.length > 1) models.forEach((model, i) => {
      const axis = model.linked_axis;
      if (!axis?.features.length) return;
      const row = append(root, "div", "", { class: "motif-form" });
      const choice = append(row, "select", "", { "aria-label": `${model.short_label} recorded feature` });
      append(choice, "option", `${model.short_label} · no recorded feature`, { value: "" });
      axis.features.forEach((item, j) => append(choice, "option",
        `${item.label} · ${item.start}–${item.end}`, { value: j }));
      for (const [kind, initial] of [["domain", "#6332ac"], ["motif", "#a92370"]]) {
        const label = append(row, "label", kind === "domain" ? "Domain colour" : "Recorded motif colour");
        const picker = append(label, "input", "", { type: "color", value: initial });
        picker.onchange = () => onLandmarkColour(i, kind, picker.value);
      }
      choice.onchange = () => {
        const item = axis.features[choice.value];
        onLandmark(i, item ? { ...item,
          toolPositions: mappedResidues(model.residues, axis.field, [item])
            .filter(position => states[i].atoms.has(position)) } : null);
      };
    });
    return { motifHits, mappedResidues };
  }

  global.pocketStructureFeatures = { mount, motifHits, mappedResidues };
})(globalThis);
