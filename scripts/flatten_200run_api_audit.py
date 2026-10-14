"""Flatten archived 200-run API response metadata without publishing response text."""

from __future__ import annotations

from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import zipfile


archive = Path(sys.argv[1])
out_dir = Path(sys.argv[2])
out_dir.mkdir(parents=True, exist_ok=True)


def utc(epoch):
    if epoch is None or epoch == "":
        return ""
    return datetime.fromtimestamp(int(epoch), timezone.utc).isoformat()


rows = []
with zipfile.ZipFile(archive) as z:
    names = sorted(n for n in z.namelist() if "/section_records/" in n and n.endswith(".json"))
    for name in names:
        record = json.loads(z.read(name))
        raw = record.get("raw_response") or {}
        usage = record.get("usage") or {}
        output = record.get("output_text") or ""
        rows.append({
            "experiment": "C2DP40_200_fixed_evidence_synthesis_runs",
            "run": record.get("run", ""),
            "section": record.get("section", ""),
            "requested_model": record.get("requested_model", ""),
            "recorded_returned_model": record.get("returned_model", ""),
            "raw_api_returned_model": raw.get("model", ""),
            "recorded_response_id": record.get("response_id", ""),
            "raw_api_response_id": raw.get("id", ""),
            "raw_api_created_epoch_seconds": raw.get("created", ""),
            "raw_api_created_utc": utc(raw.get("created")),
            "client_started_at_utc": record.get("started_at_utc", ""),
            "client_completed_at_utc": record.get("completed_at_utc", ""),
            "temperature_sent": record.get("temperature", ""),
            "thinking_sent": record.get("thinking", ""),
            "max_tokens_sent": record.get("max_tokens", ""),
            "prompt_tokens_returned": usage.get("prompt_tokens", ""),
            "completion_tokens_returned": usage.get("completion_tokens", ""),
            "total_tokens_returned": usage.get("total_tokens", ""),
            "finish_reason": record.get("finish_reason", ""),
            "prompt_sha256": record.get("prompt_sha256", ""),
            "output_sha256_computed": hashlib.sha256(output.encode("utf-8")).hexdigest(),
            "archive_record_path": name,
        })

csv_path = out_dir / "API_200run_800calls_raw_model_time_parameters.csv"
with csv_path.open("w", newline="", encoding="utf-8") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

summary = {
    "archive": archive.name,
    "record_count": len(rows),
    "run_count": len({row["run"] for row in rows}),
    "sections": dict(Counter(row["section"] for row in rows)),
    "raw_returned_models": dict(Counter(row["raw_api_returned_model"] for row in rows)),
    "all_raw_ids_unique": len({row["raw_api_response_id"] for row in rows}) == len(rows),
    "all_raw_models_match_recorded": all(row["raw_api_returned_model"] == row["recorded_returned_model"] for row in rows),
    "all_raw_ids_match_recorded": all(row["raw_api_response_id"] == row["recorded_response_id"] for row in rows),
    "raw_provider_created_present": sum(bool(row["raw_api_created_epoch_seconds"]) for row in rows),
    "client_start_present": sum(bool(row["client_started_at_utc"]) for row in rows),
    "client_completion_present": sum(bool(row["client_completed_at_utc"]) for row in rows),
    "temperature_values": sorted({str(row["temperature_sent"]) for row in rows}),
    "thinking_values": sorted({str(row["thinking_sent"]) for row in rows}),
    "max_tokens_values": sorted({str(row["max_tokens_sent"]) for row in rows}),
    "limitations": [
        "This manifest covers only the 200-run C2DP40 fixed-evidence stability experiment.",
        "Client start/completion times are not provider-returned timestamps; raw_api_created is provider-returned.",
        "Only explicitly recorded request parameters are listed; unspecified provider defaults cannot be reconstructed.",
        "Knowledge cutoff and context-window size are not in the archived response payload.",
    ],
}
(out_dir / "API_200run_800calls_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
