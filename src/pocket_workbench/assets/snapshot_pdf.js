/* Package the current rendered scene as a one-page PDF, entirely in the browser. */
((global) => {
  "use strict";
  const encoder = new TextEncoder();

  /** Escape text for a PDF literal string using a built-in font. */
  function literal(value) {
    return String(value).replace(/[\u2010-\u2015]/g, "-").replace(/·/g, " / ")
      .replace(/[^\x20-\x7e]/g, "?").replace(/[\\()]/g, "\\$&");
  }

  /** Break a short legend at a word boundary where possible. */
  function legendLines(value) {
    const words = value.replace(/·/g, "/").trim().split(/\s+/);
    const lines = [""];
    for (const word of words) {
      if (lines.at(-1).length + word.length + 1 > 36 && lines.at(-1)) lines.push("");
      lines[lines.length - 1] += `${lines.at(-1) ? " " : ""}${word}`;
    }
    if (lines.length > 2) lines[1] += "...";
    return lines.slice(0, 2).map(line => line.length > 39 ? `${line.slice(0, 36)}...` : line);
  }

  /** Encode an embedded RGB JPEG and a vector title and legend with a byte-accurate xref. */
  function encode({ jpeg, width, height, title, subtitle, legend }) {
    if (!Number.isInteger(width) || !Number.isInteger(height) || width < 1 || height < 1 ||
      !(jpeg instanceof Uint8Array) || jpeg[0] !== 255 || jpeg[1] !== 216) {
      throw new Error("The scene image is unavailable.");
    }
    const paperWidth = 842, paperHeight = 595, margin = 32;
    const imageBox = { x: 32, y: 44, width: 552, height: 464 };
    const factor = Math.min(imageBox.width / width, imageBox.height / height);
    const w = (width * factor).toFixed(2), h = (height * factor).toFixed(2);
    const x = (imageBox.x + (imageBox.width - w) / 2).toFixed(2);
    const y = (imageBox.y + (imageBox.height - h) / 2).toFixed(2);
    const commands = [];
    const text = (value, x, y, size) => commands.push(
      `0.08 0.18 0.23 rg BT /F1 ${size} Tf 1 0 0 1 ${x} ${y} Tm (${literal(value)}) Tj ET`
    );
    text(title.slice(0, 65), margin, 548, 18);
    text(subtitle.length > 119 ? `${subtitle.slice(0, 116)}...` : subtitle, margin, 527, 9);
    commands.push(`q ${w} 0 0 ${h} ${x} ${y} cm /Scene Do Q`);
    text("Scene selections", 608, 497, 12);
    if (!legend.length) text("No selections", 608, 474, 9);
    const summary = legend.length > 15 ? legend.slice(0, 14) : legend;
    for (const [index, row] of summary.entries()) {
      const y = 475 - index * 27;
      const colour = /^#[a-f0-9]{6}$/i.test(row.colour) ? row.colour : "#526b7a";
      const rgb = colour.slice(1).match(/../g).map(pair => (parseInt(pair, 16) / 255).toFixed(4));
      commands.push(`${rgb.join(" ")} rg 608 ${y - 4} 8 8 re f`);
      legendLines(String(row.text)).forEach((line, lineIndex) => text(line, 622, y - lineIndex * 11, 8.5));
    }
    if (legend.length > summary.length) text(`${legend.length - summary.length} more selections in the viewer`, 608, 62, 8);
    text("Pocket Workbench | Current view | Highlights show model evidence", margin, 21, 8);
    const contents = encoder.encode(commands.join("\n"));
    const parts = [], offsets = [0];
    let total = 0;
    const add = value => {
      const bytes = typeof value === "string" ? encoder.encode(value) : value;
      parts.push(bytes); total += bytes.byteLength;
    };
    function object(index, body, bytes) {
      offsets[index] = total;
      add(`${index} 0 obj\n${body}`);
      if (bytes) { add("\nstream\n"); add(bytes); add("\nendstream"); }
      add("\nendobj\n");
    }
    add("%PDF-1.4\n");
    add(new Uint8Array([37, 226, 227, 207, 211, 10]));
    object(1, "<< /Type /Catalog /Pages 2 0 R >>");
    object(2, "<< /Type /Pages /Kids [3 0 R] /Count 1 >>");
    object(3, `<< /Type /Page /Parent 2 0 R /MediaBox [0 0 ${paperWidth} ${paperHeight}] `
      + "/Resources << /Font << /F1 4 0 R >> /XObject << /Scene 5 0 R >> >> /Contents 6 0 R >>");
    object(4, "<< /Type /Font /Subtype /Type1 /BaseFont /Courier /Encoding /WinAnsiEncoding >>");
    object(5, `<< /Type /XObject /Subtype /Image /Width ${width} /Height ${height} `
      + `/ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /DCTDecode /Length ${jpeg.length} >>`, jpeg);
    object(6, `<< /Length ${contents.length} >>`, contents);
    const xref = total;
    add("xref\n0 7\n0000000000 65535 f \n");
    for (const offset of offsets.slice(1)) add(`${String(offset).padStart(10, "0")} 00000 n \n`);
    add(`trailer\n<< /Size 7 /Root 1 0 R >>\nstartxref\n${xref}\n%%EOF\n`);
    return new Blob(parts, { type: "application/pdf" });
  }

  /** Flatten the current snapshot onto white before converting it to JPEG for PDF. */
  async function make({ uri, title, subtitle, legend }) {
    if (!uri?.startsWith("data:image/png")) throw new Error("No rendered PNG is available.");
    const image = new Image();
    image.src = uri;
    await image.decode();
    const canvas = document.createElement("canvas");
    const scale = Math.min(1, 3000 / Math.max(image.width, image.height));
    canvas.width = Math.max(1, Math.round(image.width * scale));
    canvas.height = Math.max(1, Math.round(image.height * scale));
    const context = canvas.getContext("2d");
    context.fillStyle = "white";
    context.fillRect(0, 0, canvas.width, canvas.height);
    context.drawImage(image, 0, 0, canvas.width, canvas.height);
    const imageBlob = await new Promise(resolve => canvas.toBlob(resolve, "image/jpeg", 0.95));
    if (!imageBlob) throw new Error("Could not encode the rendered model.");
    return encode({ jpeg: new Uint8Array(await imageBlob.arrayBuffer()), width: canvas.width,
      height: canvas.height, title, subtitle, legend });
  }

  global.pocketSnapshotPDF = { make, encode };
})(globalThis);
