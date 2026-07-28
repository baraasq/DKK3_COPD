from __future__ import annotations

import csv
import gzip
import importlib.util
import io
import sys
import tarfile
import tempfile
import unittest
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from common import load_config, normalize_condition, resolve_column


def load_audit_module():
    path = ROOT / "scripts" / "00_audit_inputs.py"
    spec = importlib.util.spec_from_file_location("audit_inputs", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load audit script.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_spatial_audit_module():
    path = ROOT / "scripts" / "03_healthy_spatial_dkk3.py"
    spec = importlib.util.spec_from_file_location("healthy_spatial_dkk3", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load healthy spatial audit script.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_geomx_audit_module():
    path = ROOT / "scripts" / "03_geomx_qc" / "00_audit_gse292993_inputs.py"
    spec = importlib.util.spec_from_file_location("geomx_audit", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load GeoMx audit script.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def csv_text(fieldnames: list[str], rows: list[dict]) -> str:
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


class InputContractTests(unittest.TestCase):
    def setUp(self):
        self.config = load_config()

    def test_condition_normalization(self):
        self.assertEqual(normalize_condition("COPD", self.config), "COPD")
        self.assertEqual(normalize_condition("Healthy", self.config), "Control")
        self.assertEqual(normalize_condition("NoCLD", self.config), "Control")

    def test_case_insensitive_column_resolution(self):
        resolved = resolve_column(
            ["Donor", "Disease", "cellsubtype"],
            ["donor"],
            label="donor",
        )
        self.assertEqual(resolved, "Donor")

    def test_spatial_feature_name_decoding(self):
        spatial_module = load_spatial_audit_module()
        self.assertEqual(
            spatial_module.decode_feature_names([b"DKK3", "DKK2"]),
            ["DKK3", "DKK2"],
        )

    def test_scrinshot_audit_finds_all_sections_and_missing_dkk3(self):
        audit_module = load_audit_module()
        fields = [
            "Donor",
            "Sample",
            "Disease",
            "X",
            "Y",
            "cellclass",
            "celltype",
            "cellsubtype",
            "COL1A2",
            "TIMP1",
        ]
        with tempfile.TemporaryDirectory() as directory:
            archive_path = Path(directory) / "COPD maps and data.zip"
            with zipfile.ZipFile(
                archive_path, mode="w", compression=zipfile.ZIP_DEFLATED
            ) as archive:
                archive.writestr(
                    "COPD lung/Cell maps and histology/COPD-D1.csv",
                    csv_text(
                        fields,
                        [
                            {
                                "Donor": "D1",
                                "Sample": "S1",
                                "Disease": "COPD",
                                "X": 1,
                                "Y": 2,
                                "cellclass": "stromal",
                                "celltype": "fibroblast",
                                "cellsubtype": "fibroblast",
                                "COL1A2": 2,
                                "TIMP1": 1,
                            }
                        ],
                    ),
                )
                archive.writestr(
                    "COPD lung/Cell maps and histology/Additional maps/"
                    "Healthy-D2.csv",
                    csv_text(
                        fields,
                        [
                            {
                                "Donor": "D2",
                                "Sample": "S2",
                                "Disease": "Healthy",
                                "X": 3,
                                "Y": 4,
                                "cellclass": "epithelial",
                                "celltype": "AT2",
                                "cellsubtype": "AT2",
                                "COL1A2": 0,
                                "TIMP1": 1,
                            }
                        ],
                    ),
                )
                archive.writestr(
                    "COPD lung/Dot coordinates and DAPI/S1/"
                    "Dot-coordinates/COL1A2.csv",
                    "x,y\n1,2\n",
                )

            result, sections, genes = audit_module.audit_scrinshot(
                archive_path, self.config
            )

        self.assertEqual(result["status"], "ok")
        self.assertEqual(result["cell_map_csvs"], 2)
        self.assertEqual(len(sections), 2)
        self.assertEqual(result["sections_by_condition"]["COPD"], 1)
        self.assertEqual(result["sections_by_condition"]["Control"], 1)
        self.assertIn("COL1A2", genes)
        self.assertFalse(result["direct_spatial_gene_analysis_supported"])

    def test_geomx_audit_finds_pkc_and_dcc_dkk3(self):
        geomx_module = load_geomx_audit_module()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pkc_path = root / "panel.pkc.gz"
            with gzip.open(pkc_path, "wt", encoding="utf-8", newline="\n") as handle:
                handle.write("TargetName,CodeClass\n")
                handle.write("DKK3,Endogenous\n")
                handle.write("ACTB,Housekeeping\n")

            tar_path = root / "GSE292993_RAW.tar"
            dcc_text = "\n".join(
                [
                    "Code_Summary",
                    "TargetName,CodeClass,Count",
                    "DKK3,Endogenous,17",
                    "Negative1,Negative,2",
                    "",
                ]
            )
            with tarfile.open(tar_path, mode="w") as archive:
                payload = dcc_text.encode("utf-8")
                info = tarfile.TarInfo("nested/DSP-TEST-A01.dcc")
                info.size = len(payload)
                archive.addfile(info, io.BytesIO(payload))

            pkc = geomx_module.inspect_pkc(pkc_path, "DKK3")
            dcc, rows = geomx_module.inspect_dcc_inputs(tar_path, root / "dcc", "DKK3")

        self.assertTrue(pkc["gene_present"])
        self.assertEqual(dcc["dcc_count"], 1)
        self.assertEqual(dcc["dccs_with_gene"], 1)
        self.assertEqual(rows[0]["roi_id_guess"], "DSP-TEST-A01")
        self.assertTrue(rows[0]["contains_negative_probe_text"])


if __name__ == "__main__":
    unittest.main()
