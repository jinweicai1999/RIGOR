#!/usr/bin/env python3
import io
import os
import re
import sys
import tokenize
from pathlib import Path
from dotenv import dotenv_values
from openai import OpenAI

ROOT = Path(__file__).resolve().parents[1]

HAN_RE = re.compile(r"[\u4e00-\u9fff]")
PLACEHOLDER_RE = re.compile(r"\{[^{}]*\}")

SKIP_SUFFIX = {
    ".png", ".jpg", ".jpeg", ".gif", ".svg", ".ai", ".zip", ".pdf", ".xlsx", ".fa", ".gz", ".bed", ".csv"
}

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

class Translator:
    def __init__(self):
        cfg = dotenv_values(ROOT / ".env")
        self.client = OpenAI(api_key=cfg.get("LLM_API_KEY"), base_url=cfg.get("LLM_BASE_URL"))
        self.model = "LongCat-Flash-Lite"
        self.cache = {}

    def translate_text(self, text: str) -> str:
        if not text or not HAN_RE.search(text):
            return text
        if text in self.cache:
            return self.cache[text]

        placeholders = {}
        def _ph(m):
            key = f"__PH_{len(placeholders)}__"
            placeholders[key] = m.group(0)
            return key

        protected = PLACEHOLDER_RE.sub(_ph, text)

        prompt = (
            "Translate Chinese text to natural English. "
            "Do not change placeholders like __PH_n__. "
            "Do not add explanations. Return translation only.\n\n"
            f"TEXT:\n{protected}"
        )

        try:
            rsp = self.client.chat.completions.create(
                model=self.model,
                messages=[{"role": "user", "content": prompt}],
                temperature=0,
            )
            out = (rsp.choices[0].message.content or "").strip()
            if not out:
                out = text
        except Exception:
            out = text

        for k, v in placeholders.items():
            out = out.replace(k, v)

        self.cache[text] = out
        return out


def contains_han(s: str) -> bool:
    return bool(HAN_RE.search(s or ""))


def process_markdown_or_text(path: Path, tr: Translator) -> bool:
    text = path.read_text(encoding="utf-8", errors="ignore")
    lines = text.splitlines(keepends=True)
    changed = False
    out_lines = []
    for line in lines:
        if contains_han(line):
            nl = tr.translate_text(line.rstrip("\n"))
            out_lines.append(nl + ("\n" if line.endswith("\n") else ""))
            changed = changed or (nl != line.rstrip("\n"))
        else:
            out_lines.append(line)
    if changed:
        path.write_text("".join(out_lines), encoding="utf-8")
    return changed


def split_string_token(tok: str):
    m = re.match(r"(?is)^([rubf]*)('''|\"\"\"|'|\")(.*)(\2)$", tok)
    if not m:
        return None
    return m.group(1), m.group(2), m.group(3), m.group(4)


def process_python(path: Path, tr: Translator) -> bool:
    text = path.read_text(encoding="utf-8", errors="ignore")
    stream = io.StringIO(text)
    toks = list(tokenize.generate_tokens(stream.readline))
    changed = False
    new_toks = []

    for tok in toks:
        ttype, tstr, start, end, line = tok
        new_str = tstr

        if ttype == tokenize.COMMENT and contains_han(tstr):
            body = tstr[1:]
            body_t = tr.translate_text(body)
            new_str = "#" + body_t

        elif ttype == tokenize.STRING and contains_han(tstr):
            parsed = split_string_token(tstr)
            if parsed:
                prefix, q1, body, q2 = parsed
                body_t = tr.translate_text(body)
                new_str = f"{prefix}{q1}{body_t}{q2}"

        if new_str != tstr:
            changed = True
            new_toks.append((ttype, new_str))
        else:
            new_toks.append((ttype, tstr))

    if changed:
        rebuilt = tokenize.untokenize(new_toks)
        path.write_text(rebuilt, encoding="utf-8")
    return changed


def main():
    tr = Translator()
    changed_files = []

    for rel in TARGET_FILES:
        path = ROOT / rel
        if not path.exists() or path.suffix.lower() in SKIP_SUFFIX:
            continue
        try:
            if path.suffix.lower() == ".py":
                changed = process_python(path, tr)
            else:
                changed = process_markdown_or_text(path, tr)
            if changed:
                changed_files.append(rel)
                print(f"UPDATED {rel}")
            else:
                print(f"SKIP {rel}")
        except Exception as e:
            print(f"ERROR {rel}: {e}")

    print("\n=== Done ===")
    print(f"Changed files: {len(changed_files)}")
    for f in changed_files:
        print(f)


if __name__ == "__main__":
    main()
