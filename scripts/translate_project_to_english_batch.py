#!/usr/bin/env python3
import io
import json
import re
import tokenize
from pathlib import Path
from collections import OrderedDict
from dotenv import dotenv_values
from openai import OpenAI

ROOT = Path(__file__).resolve().parents[1]
HAN_RE = re.compile(r"[\u4e00-\u9fff]")

TARGET_FILES = [
    "PROJECT_SUMMARY.md",
    "TRANSLATION_SUMMARY.md",
    "apps/lncrna_pipeline.py",
    "apps/mirna_pipeline.py",
    "apps/protein_batch_pipeline.py",
    "apps/protein_function_pipeline.py",
    "apps/protein_trait_pipeline_legacy.py",
    "apps/sequence_identifier.py",
    "apps/web_dashboard.py",
    "docs/MiRNA_Analysis.md",
    "docs/Protein_Coding_Gene.md",
    "docs/lncrna_pipeline_updated_logic.md",
    "docs/mirna_sequence_mode.md",
    "docs/project_logic.md",
    "docs/project_overview.md",
    "scripts/build_mirna_target_db.py",
    "src/agents/gwas_agent.py",
    "src/agents/literature_agent.py",
    "src/agents/lncrna_agent.py",
    "src/agents/macro_network_agent.py",
    "src/agents/mirna_sequence_agent.py",
    "src/agents/ncrna_agent.py",
    "src/utils/coexpression_engine.py",
    "src/utils/cross_omics_engine.py",
    "src/utils/epigenomic_analyzer.py",
    "src/utils/expression_analyzer.py",
    "src/utils/expression_pattern_service.py",
    "src/utils/gene_annotator.py",
    "src/utils/genome_consistency_hint_service.py",
    "src/utils/gff_explorer.py",
    "src/utils/known_lncrna_annotator.py",
    "src/utils/lncrna_annotator.py",
    "src/utils/parse_plink.py",
    "src/utils/parse_real_gwas.py",
    "src/utils/protein_batch_analysis.py",
    "src/utils/protein_feature_service.py",
    "src/utils/rna_type_detector.py",
    "translate_chinese.py",
]


def has_han(s: str) -> bool:
    return bool(HAN_RE.search(s or ""))


def split_string_token(tok: str):
    m = re.match(r"(?is)^([rubf]*)('''|\"\"\"|'|\")(.*)(\2)$", tok)
    if not m:
        return None
    return m.group(1), m.group(2), m.group(3), m.group(4)


def collect_segments(path: Path):
    text = path.read_text(encoding="utf-8", errors="ignore")
    segs = []
    if path.suffix.lower() == ".py":
        toks = list(tokenize.generate_tokens(io.StringIO(text).readline))
        for ttype, tstr, *_ in toks:
            if ttype == tokenize.COMMENT and has_han(tstr):
                segs.append(("comment", tstr[1:]))
            elif ttype == tokenize.STRING and has_han(tstr):
                parsed = split_string_token(tstr)
                if parsed:
                    _, _, body, _ = parsed
                    segs.append(("string", body))
    else:
        for line in text.splitlines():
            if has_han(line):
                segs.append(("line", line))
    return segs


def batch_translate(unique_segments):
    cfg = dotenv_values(ROOT / ".env")
    client = OpenAI(api_key=cfg.get("LLM_API_KEY"), base_url=cfg.get("LLM_BASE_URL"), timeout=60.0)
    out = {}

    items = list(unique_segments)
    chunk_size = 30
    for i in range(0, len(items), chunk_size):
        chunk = items[i:i+chunk_size]
        payload = {str(j): s for j, s in enumerate(chunk)}
        prompt = (
            "Translate each Chinese text value to natural English. Keep placeholders like {name} unchanged. "
            "Return strict JSON only with the same keys.\n\n"
            + json.dumps(payload, ensure_ascii=False)
        )
        try:
            rsp = client.chat.completions.create(
                model="LongCat-Flash-Lite",
                messages=[{"role": "user", "content": prompt}],
                temperature=0,
            )
            txt = (rsp.choices[0].message.content or "").strip()
            m = re.search(r"\{[\s\S]*\}", txt)
            if not m:
                for s in chunk:
                    out[s] = s
                continue
            data = json.loads(m.group(0))
            for k, src in payload.items():
                out[src] = data.get(k, src)
        except Exception:
            for s in chunk:
                out[s] = s
        print(f"translated batch {i//chunk_size + 1}/{(len(items)+chunk_size-1)//chunk_size}", flush=True)
    return out


def apply_file(path: Path, tr_map):
    text = path.read_text(encoding="utf-8", errors="ignore")
    changed = False

    if path.suffix.lower() == ".py":
        toks = list(tokenize.generate_tokens(io.StringIO(text).readline))
        new_toks = []
        for ttype, tstr, *_ in toks:
            new_str = tstr
            if ttype == tokenize.COMMENT and has_han(tstr):
                body = tstr[1:]
                body_t = tr_map.get(body, body)
                new_str = "#" + body_t
            elif ttype == tokenize.STRING and has_han(tstr):
                parsed = split_string_token(tstr)
                if parsed:
                    p, q1, body, q2 = parsed
                    body_t = tr_map.get(body, body)
                    new_str = f"{p}{q1}{body_t}{q2}"
            if new_str != tstr:
                changed = True
            new_toks.append((ttype, new_str))
        if changed:
            text = tokenize.untokenize(new_toks)
            path.write_text(text, encoding="utf-8")
        return changed

    lines = text.splitlines(keepends=True)
    out_lines = []
    for line in lines:
        body = line[:-1] if line.endswith("\n") else line
        if has_han(body):
            new_body = tr_map.get(body, body)
            if new_body != body:
                changed = True
            out_lines.append(new_body + ("\n" if line.endswith("\n") else ""))
        else:
            out_lines.append(line)
    if changed:
        path.write_text("".join(out_lines), encoding="utf-8")
    return changed


def main():
    segment_set = OrderedDict()
    files = []
    for rel in TARGET_FILES:
        p = ROOT / rel
        if not p.exists():
            continue
        files.append(p)
        for _, seg in collect_segments(p):
            if has_han(seg):
                segment_set.setdefault(seg, None)

    uniq = list(segment_set.keys())
    print(f"unique segments: {len(uniq)}", flush=True)
    if not uniq:
        print("nothing to translate", flush=True)
        return

    tr_map = batch_translate(uniq)

    changed_count = 0
    for p in files:
        try:
            if apply_file(p, tr_map):
                changed_count += 1
                print(f"UPDATED {p.relative_to(ROOT)}", flush=True)
        except Exception as e:
            print(f"ERROR {p.relative_to(ROOT)}: {e}", flush=True)

    print(f"done changed_files={changed_count}", flush=True)


if __name__ == "__main__":
    main()
