#!/usr/bin/env python3
"""Run RIGOR single-gene prediction and assemble a Markdown/PDF report."""

from __future__ import annotations

import argparse
import base64
import json
import os
import re
import subprocess
import sys
import textwrap
from datetime import datetime
from pathlib import Path
from typing import Any

from PIL import Image, ImageDraw, ImageFont


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PIPELINE = PROJECT_ROOT / "apps" / "protein_function_pipeline.py"
STEP_BEGIN = "PROTEIN_STEP_JSON_BEGIN"
STEP_END = "PROTEIN_STEP_JSON_END"
RESULT_BEGIN = "PROTEIN_RESULT_JSON_BEGIN"
RESULT_END = "PROTEIN_RESULT_JSON_END"


def extract_json_blocks(text: str, begin: str, end: str) -> list[dict[str, Any]]:
    pattern = re.compile(
        rf"{re.escape(begin)}\s*(.*?)\s*{re.escape(end)}",
        re.DOTALL,
    )
    blocks: list[dict[str, Any]] = []
    for match in pattern.finditer(text):
        raw = match.group(1).strip()
        try:
            blocks.append(json.loads(raw))
        except json.JSONDecodeError:
            continue
    return blocks


def save_b64_png(value: str, path: Path) -> str:
    if not value:
        return ""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(base64.b64decode(value))
    return path.name


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_")


def run_gene(gene: str, out_dir: Path, env: dict[str, str], force: bool = False) -> dict[str, Any]:
    raw_dir = out_dir / "raw"
    stdout_path = raw_dir / f"{safe_name(gene)}.stdout.log"
    stderr_path = raw_dir / f"{safe_name(gene)}.stderr.log"
    if stdout_path.exists() and not force:
        stdout = stdout_path.read_text(encoding="utf-8", errors="replace")
        stderr = stderr_path.read_text(encoding="utf-8", errors="replace") if stderr_path.exists() else ""
        steps = extract_json_blocks(stdout, STEP_BEGIN, STEP_END)
        results = extract_json_blocks(stdout, RESULT_BEGIN, RESULT_END)
        return {
            "gene": gene,
            "returncode": 0 if results else 1,
            "steps": steps,
            "result": results[-1] if results else {},
            "stderr": stderr,
        }

    cmd = [sys.executable, str(PIPELINE), "--gene", gene]
    proc = subprocess.run(
        cmd,
        cwd=PROJECT_ROOT,
        env=env,
        text=True,
        capture_output=True,
        timeout=1800,
    )
    raw_dir.mkdir(parents=True, exist_ok=True)
    stdout_path.write_text(proc.stdout, encoding="utf-8")
    stderr_path.write_text(proc.stderr, encoding="utf-8")

    steps = extract_json_blocks(proc.stdout, STEP_BEGIN, STEP_END)
    results = extract_json_blocks(proc.stdout, RESULT_BEGIN, RESULT_END)
    return {
        "gene": gene,
        "returncode": proc.returncode,
        "steps": steps,
        "result": results[-1] if results else {},
        "stderr": proc.stderr,
    }


def step_by_name(steps: list[dict[str, Any]], name: str) -> dict[str, Any]:
    for step in steps:
        if step.get("step") == name:
            return step.get("data", {})
    return {}


def first_lines(text: str, max_lines: int = 18) -> str:
    lines = [line.rstrip() for line in (text or "").splitlines() if line.strip()]
    return "\n".join(lines[:max_lines])


def iter_sections(sections: Any) -> list[tuple[str, str]]:
    def clean_pair(title: Any, content: Any) -> tuple[str, str] | None:
        title_s = str(title)
        content_s = str(content or "").strip()
        if content_s.startswith("{") and content_s.endswith("}"):
            try:
                parsed = json.loads(content_s)
                if isinstance(parsed, dict):
                    title_s = str(parsed.get("title") or title_s)
                    content_s = str(parsed.get("content") or parsed.get("text") or parsed.get("body") or "").strip()
            except json.JSONDecodeError:
                pass
        title_s = "Integrated evidence summary" if title_s == "Untitled Section" else title_s
        if not content_s:
            return None
        return title_s, content_s

    if isinstance(sections, dict):
        pairs = []
        for key, value in sections.items():
            pair = clean_pair(key, value)
            if pair:
                pairs.append(pair)
        return pairs
    if isinstance(sections, list):
        normalized = []
        for idx, item in enumerate(sections, start=1):
            if isinstance(item, str):
                stripped = item.strip()
                if stripped.startswith("{") and stripped.endswith("}"):
                    try:
                        item = json.loads(stripped)
                    except json.JSONDecodeError:
                        pass
            if isinstance(item, dict):
                title = item.get("title") or item.get("heading") or f"Section {idx}"
                content = (
                    item.get("content")
                    or item.get("text")
                    or item.get("body")
                    or item.get("report")
                    or ""
                )
                if not content and len(item) > 1:
                    content = json.dumps(item, ensure_ascii=False)
            else:
                title = f"Section {idx}"
                content = str(item)
            pair = clean_pair(title, content)
            if pair:
                normalized.append(pair)
        return normalized
    if sections:
        return [("Section", str(sections))]
    return []


def save_gene_images(gene_result: dict[str, Any], image_dir: Path) -> dict[str, list[str] | str]:
    gene = gene_result["gene"]
    prefix = safe_name(gene)
    steps = gene_result["steps"]
    images: dict[str, list[str] | str] = {"expression": []}

    basic = step_by_name(steps, "basic_info")
    protein_png = save_b64_png(
        basic.get("plot_png_b64", ""),
        image_dir / f"{prefix}_protein_features.png",
    )
    if protein_png:
        images["protein"] = protein_png

    expression = step_by_name(steps, "expression")
    for idx, tr in enumerate(expression.get("transcript_results", []) or [], start=1):
        tid = safe_name(tr.get("transcript_id") or f"transcript_{idx}")
        dev = save_b64_png(
            tr.get("developmental_plot_b64", ""),
            image_dir / f"{prefix}_{tid}_developmental_expression.png",
        )
        treat = save_b64_png(
            tr.get("treatment_plot_b64", ""),
            image_dir / f"{prefix}_{tid}_treatment_expression.png",
        )
        if dev:
            images["expression"].append(dev)  # type: ignore[union-attr]
        if treat:
            images["expression"].append(treat)  # type: ignore[union-attr]

    return images


def make_markdown(all_results: list[dict[str, Any]], out_dir: Path) -> str:
    image_dir = out_dir / "images"
    date = datetime.now().strftime("%Y-%m-%d %H:%M")
    lines: list[str] = [
        "# RIGOR single-gene prediction report",
        "",
        f"- Generated: {date}",
        "- Workflow: RIGOR single-gene protein-coding inference",
        "- LLM backend: LongCat-Flash-Lite through an OpenAI-compatible API endpoint",
        "- Evidence layers: annotation/domain evidence, homologous literature, co-expression neighborhood, companion-gene literature, promoter epigenomic state, MSU expression pattern, and GWAS trait context.",
        "",
        "> Note: API credentials are intentionally omitted from this report.",
        "",
    ]

    for item in all_results:
        gene = item["gene"]
        steps = item["steps"]
        result = item["result"]
        images = save_gene_images(item, image_dir)
        basic = step_by_name(steps, "basic_info")
        homology = step_by_name(steps, "homology")
        coexp = step_by_name(steps, "coexpression")
        epi = step_by_name(steps, "epigenomic")
        expr = step_by_name(steps, "expression")
        gwas = step_by_name(steps, "gwas")
        final_sections = result.get("final_sections", {})
        gwas_sections = result.get("gwas_sections", {})

        lines.extend(
            [
                f"## {gene}",
                "",
                "### Basic annotation",
                "",
                f"- RAP ID: {basic.get('target') or result.get('target') or gene}",
                f"- Symbol: {basic.get('symbol') or result.get('symbol') or 'N/A'}",
                f"- Description: {basic.get('description') or result.get('description') or 'N/A'}",
                f"- Genomic location: {basic.get('location') or result.get('location') or 'N/A'}",
                f"- UniProt: {basic.get('uniprot') or 'N/A'}",
                f"- Protein length: {basic.get('protein_length') or 'N/A'}",
                "",
            ]
        )
        if images.get("protein"):
            lines.extend(
                [
                    f"![Protein features](images/{images['protein']})",
                    "",
                ]
            )

        interpro = basic.get("interpro_domains") or []
        pfam = basic.get("pfam_domains") or []
        if interpro or pfam:
            lines.extend(["### Domain and protein-feature evidence", ""])
            if interpro:
                lines.append(f"- InterPro domain fragments detected: {len(interpro)}")
            if pfam:
                lines.append(f"- Pfam domain fragments detected: {len(pfam)}")
            locs = basic.get("subcellular_locations") or []
            if locs:
                lines.append(f"- Predicted/annotated subcellular localization: {', '.join(locs)}")
            lines.append("")

        lines.extend(["### Homology-guided literature inference", ""])
        lines.append(first_lines(homology.get("report", ""), 20) or "No homology literature inference was returned.")
        refs = homology.get("references") or []
        if refs:
            lines.append("")
            lines.append("Key retrieved references:")
            for ref in refs[:8]:
                title = ref.get("title") or "Untitled"
                pmid = ref.get("pmid") or "N/A"
                lines.append(f"- PMID {pmid}: {title}")
        lines.append("")

        lines.extend(["### Co-expression and companion-gene evidence", ""])
        lines.append(first_lines(coexp.get("report", ""), 20) or "No co-expression inference was returned.")
        lines.append("")
        lines.append(f"- Co-expressed companion genes used: {coexp.get('network_gene_count', 0)}")
        go_terms = coexp.get("go_terms") or []
        if go_terms:
            lines.append("- Top local GO themes: " + "; ".join(f"{x.get('term')} ({x.get('count')})" for x in go_terms[:10]))
        lines.append("")

        lines.extend(["### Promoter epigenomic context", ""])
        lines.append(first_lines(epi.get("report", ""), 16) or "No promoter epigenomic inference was returned.")
        tissue_scores = epi.get("tissue_scores") or []
        if tissue_scores:
            compact = "; ".join(f"{x.get('tissue')}: {x.get('state')} ({x.get('score')})" for x in tissue_scores[:10])
            lines.append("")
            lines.append(f"- Promoter-state summary: {compact}")
        lines.append("")

        lines.extend(["### MSU expression-pattern context", ""])
        lines.append(expr.get("summary") or "No MSU expression pattern summary was returned.")
        expr_images = images.get("expression") or []
        for img in expr_images:
            lines.extend(["", f"![Expression pattern](images/{img})"])
        lines.append("")

        lines.extend(["### GWAS and breeding-oriented trait context", ""])
        if gwas_sections:
            for title, content in iter_sections(gwas_sections):
                lines.append(f"**{title}**")
                lines.append("")
                lines.append(first_lines(content, 16))
                lines.append("")
        elif gwas.get("report_sections"):
            for title, content in iter_sections(gwas.get("report_sections")):
                lines.append(f"**{title}**")
                lines.append("")
                lines.append(first_lines(content, 16))
                lines.append("")
        else:
            lines.append("No GWAS summary was returned.")
            lines.append("")

        lines.extend(["### Integrated prediction", ""])
        if final_sections:
            for title, content in iter_sections(final_sections):
                lines.append(f"**{title}**")
                lines.append("")
                lines.append(first_lines(content, 18))
                lines.append("")
        else:
            lines.append("No integrated final section was returned.")
            lines.append("")

        if item["returncode"] != 0:
            lines.extend(
                [
                    "### Runtime warning",
                    "",
                    f"- The pipeline exited with return code {item['returncode']}. See raw logs for details.",
                    "",
                ]
            )

    return "\n".join(lines).rstrip() + "\n"


def load_font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/Library/Fonts/Arial Unicode.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ]
    for candidate in candidates:
        if Path(candidate).exists():
            return ImageFont.truetype(candidate, size=size)
    return ImageFont.load_default()


def markdown_to_pdf(markdown_text: str, pdf_path: Path) -> None:
    """Create a readable PDF from Markdown using PIL only."""
    page_w, page_h = 1654, 2339
    margin = 105
    line_gap = 12
    font = load_font(30)
    h1 = load_font(48)
    h2 = load_font(40)
    h3 = load_font(34)
    small = load_font(24)
    pages: list[Image.Image] = []
    page = Image.new("RGB", (page_w, page_h), "white")
    draw = ImageDraw.Draw(page)
    y = margin

    def new_page() -> None:
        nonlocal page, draw, y
        pages.append(page)
        page = Image.new("RGB", (page_w, page_h), "white")
        draw = ImageDraw.Draw(page)
        y = margin

    def add_text(raw: str, selected_font: ImageFont.ImageFont, fill: str = "#1f2933") -> None:
        nonlocal y
        if not raw.strip():
            y += 22
            return
        width_chars = max(40, int((page_w - 2 * margin) / (selected_font.size * 0.55 if hasattr(selected_font, "size") else 16)))
        prefix = ""
        text = raw
        if raw.startswith("- "):
            prefix = "- "
            text = raw[2:]
        wrapped = textwrap.wrap(text, width=width_chars, break_long_words=False) or [""]
        for idx, part in enumerate(wrapped):
            line = (prefix if idx == 0 else "  ") + part
            bbox = draw.textbbox((margin, y), line, font=selected_font)
            if y + (bbox[3] - bbox[1]) + margin > page_h:
                new_page()
            draw.text((margin, y), line, font=selected_font, fill=fill)
            y += (bbox[3] - bbox[1]) + line_gap

    for raw_line in markdown_text.splitlines():
        line = raw_line.rstrip()
        if line.startswith("!["):
            match = re.search(r"\(([^)]+)\)", line)
            if not match:
                continue
            image_path = pdf_path.parent / match.group(1)
            if not image_path.exists():
                continue
            with Image.open(image_path) as src:
                src = src.convert("RGB")
                max_w = page_w - 2 * margin
                max_h = 760
                scale = min(max_w / src.width, max_h / src.height, 1.0)
                resized = src.resize((int(src.width * scale), int(src.height * scale)))
                if y + resized.height + margin > page_h:
                    new_page()
                page.paste(resized, (margin, y))
                y += resized.height + 28
            continue
        if line.startswith("# "):
            add_text(line[2:], h1, "#0f172a")
        elif line.startswith("## "):
            y += 15
            add_text(line[3:], h2, "#123c69")
        elif line.startswith("### "):
            y += 10
            add_text(line[4:], h3, "#1b4d3e")
        elif line.startswith("> "):
            add_text(line[2:], small, "#556070")
        elif line.startswith("**") and line.endswith("**"):
            add_text(line.strip("*"), h3, "#334155")
        else:
            add_text(line, font)
    pages.append(page)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    pages[0].save(pdf_path, save_all=True, append_images=pages[1:], resolution=150.0)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genes", nargs="+", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--base-url", required=True)
    parser.add_argument("--api-key", required=True)
    parser.add_argument("--model", default="LongCat-Flash-Lite")
    parser.add_argument("--force", action="store_true", help="rerun predictions even if raw logs exist")
    args = parser.parse_args()

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env.update(
        {
            "APP_LANG": "en",
            "LLM_BASE_URL": args.base_url,
            "LLM_API_KEY": args.api_key,
            "LLM_MODEL_NAME": args.model,
            "MPLCONFIGDIR": str(out_dir / ".mplconfig"),
        }
    )
    (out_dir / ".mplconfig").mkdir(parents=True, exist_ok=True)

    results = []
    for gene in args.genes:
        print(f"[RIGOR] running {gene}", flush=True)
        results.append(run_gene(gene, out_dir, env, force=args.force))

    md = make_markdown(results, out_dir)
    md_path = out_dir / "RIGOR_single_gene_predictions.md"
    md_path.write_text(md, encoding="utf-8")
    pdf_path = out_dir / "RIGOR_single_gene_predictions.pdf"
    markdown_to_pdf(md, pdf_path)

    summary = {
        "markdown": str(md_path),
        "pdf": str(pdf_path),
        "out_dir": str(out_dir),
        "genes": [
            {
                "gene": item["gene"],
                "returncode": item["returncode"],
                "step_count": len(item["steps"]),
                "has_result": bool(item["result"]),
            }
            for item in results
        ],
    }
    (out_dir / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)
    return 0 if all(item["returncode"] == 0 and item["result"] for item in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
