#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config, project_path


DEFAULT_PREPROCESSING_CODE = (
    "intermediate/gse302339_scanpy_workflow_code/1_preprocessing_doublet_detection.py"
)

HARMONY_WRAPPER_CALL = "sce.pp.harmony_integrate(adata, 'batch', **kwargs)"

MANUAL_HARMONY_BLOCK = """# CODEx-PATCH: run HarmonyPy directly and validate orientation.
import harmonypy as hm
harmony_out = hm.run_harmony(
    adata.obsm['X_pca'],
    adata.obs,
    vars_use=['batch'],
    **kwargs,
)
z_corr = harmony_out.Z_corr
if hasattr(z_corr, 'to_numpy'):
    z_corr = z_corr.to_numpy()
z_corr = np.asarray(z_corr)
expected_shape = adata.obsm['X_pca'].shape
if z_corr.shape == expected_shape:
    adata.obsm['X_pca_harmony'] = z_corr
elif z_corr.T.shape == expected_shape:
    adata.obsm['X_pca_harmony'] = z_corr.T
else:
    raise ValueError(
        'Harmony output shape mismatch: '
        f'Z_corr={z_corr.shape}, Z_corr.T={z_corr.T.shape}, '
        f'expected={expected_shape}'
    )"""


def patch_preprocessing_code(path: Path) -> dict:
    if not path.exists():
        return {
            "path": str(path),
            "exists": False,
            "patched": False,
            "reason": "file_missing",
        }

    text = path.read_text(encoding="utf-8")
    if "CODEx-PATCH: run HarmonyPy directly" in text:
        return {
            "path": str(path),
            "exists": True,
            "patched": False,
            "reason": "already_patched",
        }
    if HARMONY_WRAPPER_CALL not in text:
        return {
            "path": str(path),
            "exists": True,
            "patched": False,
            "reason": "harmony_wrapper_call_not_found",
        }

    patched = text.replace(HARMONY_WRAPPER_CALL, MANUAL_HARMONY_BLOCK, 1)
    path.write_text(patched, encoding="utf-8")
    return {
        "path": str(path),
        "exists": True,
        "patched": True,
        "reason": "replaced_scanpy_harmony_wrapper",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Patch the exported GSE302339 preprocessing notebook code for runtime "
            "compatibility with the active Scanpy/HarmonyPy stack."
        )
    )
    parser.add_argument("--preprocessing-code", default=DEFAULT_PREPROCESSING_CODE)
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    results_dirs = ensure_results_dirs(load_config())
    path = project_path(args.preprocessing_code)

    summary = patch_preprocessing_code(path)
    summary_path = results_dirs["meta"] / "gse302339_author_preprocessing_patch_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print()
    print(summary_path)

    if args.strict and not (summary["patched"] or summary["reason"] == "already_patched"):
        raise SystemExit(
            "Strict preprocessing runtime patch failed: " + str(summary["reason"])
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
