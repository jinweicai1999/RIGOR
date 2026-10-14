"""Export literal runtime prompt expressions for source-review reproducibility.

This is a documentation extractor. The source expressions in apps/ and src/ remain
the executable source of truth; generated Markdown is not imported at runtime.
"""

from __future__ import annotations

import ast
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIRS = (ROOT / "apps", ROOT / "src")
OUTPUT = ROOT / "prompts" / "runtime_prompt_templates.md"


def assigned_names(node: ast.Assign | ast.AnnAssign) -> list[str]:
    targets = node.targets if isinstance(node, ast.Assign) else [node.target]
    return [target.id for target in targets if isinstance(target, ast.Name)]


def extract() -> list[tuple[str, int, str, str]]:
    records = []
    for source_dir in SOURCE_DIRS:
        for path in sorted(source_dir.rglob("*.py")):
            text = path.read_text(encoding="utf-8")
            tree = ast.parse(text, filename=str(path))
            for node in ast.walk(tree):
                if not isinstance(node, (ast.Assign, ast.AnnAssign)):
                    continue
                names = assigned_names(node)
                if not any(name in {"prompt", "orphan_prompt"} for name in names):
                    continue
                value = node.value
                if not isinstance(value, (ast.JoinedStr, ast.Constant)):
                    continue
                if isinstance(value, ast.Constant) and not isinstance(value.value, str):
                    continue
                expression = ast.get_source_segment(text, value)
                if expression:
                    records.append((str(path.relative_to(ROOT)), node.lineno, names[0], expression))
    return sorted(records, key=lambda row: (row[0], row[1]))


def main() -> None:
    records = extract()
    lines = [
        "# RIGOR runtime prompt expressions",
        "",
        "This file is a line-anchored export of the literal prompt expressions in the",
        "released Python source. The executable source of truth is the cited code;",
        "the formerly empty `src/prompts/system_prompts.py` is not imported by the",
        "current runtime. Dynamic evidence fields are shown as Python f-string",
        "expressions and are populated only during a run. This export contains no",
        "patient, user, API-key, or dataset content.",
        "",
        f"Extracted prompt expressions: {len(records)}.",
        "",
    ]
    for path, line, name, expression in records:
        lines.extend([f"## `{path}:{line}` — `{name}`", "", "```python", expression, "```", ""])
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text("\n".join(lines), encoding="utf-8")
    print(f"{len(records)} prompt expressions -> {OUTPUT}")


if __name__ == "__main__":
    main()
