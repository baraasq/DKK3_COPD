#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config, project_path


DEFAULT_ANNOTATION_CODE = (
    "intermediate/gse302339_scanpy_workflow_code/2_celltype_annotation.py"
)
OPTIONAL_MERGE_ARTIFACT = "output/adata_mergedmeta_abt_cr8"
PATCH_MARKER = "CODEx-PATCH: stop before optional missing ABT/meta-merge section"

OPEN_OPTIONAL_ARTIFACT_PATTERN = re.compile(
    r"""(?P<indent>^[ \t]*)with\s+open\(\s*["']output/adata_mergedmeta_abt_cr8["']\s*,\s*["']rb["']\s*\)\s+as\s+file\s*:""",
    re.MULTILINE,
)


def patch_annotation_code(path: Path) -> dict:
    if not path.exists():
        return {
            "path": str(path),
            "exists": False,
            "patched": False,
            "reason": "file_missing",
        }

    text = path.read_text(encoding="utf-8")
    if PATCH_MARKER in text:
        return {
            "path": str(path),
            "exists": True,
            "patched": False,
            "reason": "already_patched",
        }

    match = OPEN_OPTIONAL_ARTIFACT_PATTERN.search(text)
    if not match:
        return {
            "path": str(path),
            "exists": True,
            "patched": False,
            "reason": "optional_artifact_read_not_found",
        }

    indent = match.group("indent")
    patch = "\n".join(
        [
            f"{indent}# {PATCH_MARKER}",
            f"{indent}print('CODEx-PATCH: stopping before optional ABT/meta-merge section; main cell-type annotation outputs should already be available.')",
            f"{indent}raise SystemExit(0)",
            "",
        ]
    )
    patched = text[: match.start()] + patch + text[match.start() :]
    path.write_text(patched, encoding="utf-8")
    return {
        "path": str(path),
        "exists": True,
        "patched": True,
        "reason": "inserted_clean_stop_before_optional_artifact_read",
        "optional_artifact": OPTIONAL_MERGE_ARTIFACT,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Patch the exported GSE302339 cell-type annotation code to stop "
            "cleanly before an optional ABT/meta-merge section that reads an "
            "artifact not needed for parenchymal signature building."
        )
    )
    parser.add_argument("--annotation-code", default=DEFAULT_ANNOTATION_CODE)
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    results_dirs = ensure_results_dirs(load_config())
    path = project_path(args.annotation_code)

    summary = patch_annotation_code(path)
    summary_path = (
        results_dirs["meta"] / "gse302339_celltype_annotation_signature_patch_summary.json"
    )
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print()
    print(summary_path)

    if args.strict and not (summary["patched"] or summary["reason"] == "already_patched"):
        raise SystemExit(
            "Strict celltype annotation patch failed: " + str(summary["reason"])
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
