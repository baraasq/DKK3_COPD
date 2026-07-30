from __future__ import annotations

import csv
import py_compile
import gzip
import importlib.util
import io
import json
import random
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


def load_geomx_profile_module():
    path = ROOT / "scripts" / "03_geomx_qc" / "01_profile_gse292993_objects.py"
    spec = importlib.util.spec_from_file_location("geomx_profile", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load GeoMx profile script.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_geomx_dcc_qc_module():
    path = ROOT / "scripts" / "03_geomx_qc" / "02_extract_gse292993_dcc_qc.py"
    spec = importlib.util.spec_from_file_location("geomx_dcc_qc", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load GeoMx DCC QC script.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_geomx_roi_qc_module():
    path = ROOT / "scripts" / "03_geomx_qc" / "03_flag_gse292993_roi_qc.py"
    spec = importlib.util.spec_from_file_location("geomx_roi_qc", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load GeoMx ROI QC script.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_geomx_loq_module():
    path = ROOT / "scripts" / "03_geomx_qc" / "04_compute_gse292993_dkk3_loq.py"
    spec = importlib.util.spec_from_file_location("geomx_loq", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load GeoMx LOQ script.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_geomx_qc_plot_module():
    path = ROOT / "scripts" / "03_geomx_qc" / "05_plot_gse292993_roi_qc.py"
    spec = importlib.util.spec_from_file_location("geomx_qc_plot", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load GeoMx QC plot script.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_geomx_matrix_module():
    path = ROOT / "scripts" / "04_normalization" / "00_build_gse292993_geomx_wta_matrix.py"
    spec = importlib.util.spec_from_file_location("geomx_matrix", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load GeoMx WTA matrix script.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_scrna_audit_module():
    path = ROOT / "scripts" / "06_celltype" / "00_audit_scrna_reference.py"
    spec = importlib.util.spec_from_file_location("scrna_audit", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load scRNA reference audit script.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_deconv_module():
    path = ROOT / "scripts" / "06_celltype" / "01_deconvolve_gse292993_compartment_nnls.py"
    spec = importlib.util.spec_from_file_location("compartment_deconv", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load compartment deconvolution script.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_scanpy_notebook_inspection_module():
    path = ROOT / "scripts" / "06_celltype" / "02_inspect_gse302339_scanpy_notebooks.py"
    spec = importlib.util.spec_from_file_location("scanpy_notebook_inspection", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load Scanpy notebook inspection script.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_annotation_code_export_module():
    path = ROOT / "scripts" / "06_celltype" / "03_extract_gse302339_annotation_code.py"
    spec = importlib.util.spec_from_file_location("annotation_code_export", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load annotation code export script.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_author_celltype_map_module():
    path = ROOT / "scripts" / "06_celltype" / "04_extract_gse302339_author_celltype_maps.py"
    spec = importlib.util.spec_from_file_location("author_celltype_maps", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load author celltype map extraction script.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_author_signature_module():
    path = ROOT / "scripts" / "06_celltype" / "05_build_gse302339_author_signatures.py"
    spec = importlib.util.spec_from_file_location("author_signatures", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load author signature builder script.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_author_input_helper_module():
    path = ROOT / "scripts" / "06_celltype" / "06_prepare_gse302339_author_input_helpers.py"
    spec = importlib.util.spec_from_file_location("author_input_helpers", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load author input helper script.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_dkk3_summary_module():
    path = ROOT / "scripts" / "05_dkk3" / "00_summarize_gse292993_dkk3.py"
    spec = importlib.util.spec_from_file_location("dkk3_summary", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load DKK3 summary script.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_dkk3_effect_module():
    path = ROOT / "scripts" / "05_dkk3" / "01_test_gse292993_dkk3_donor_effects.py"
    spec = importlib.util.spec_from_file_location("dkk3_effects", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load DKK3 effect script.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_dkk3_strata_module():
    path = ROOT / "scripts" / "05_dkk3" / "02_test_gse292993_dkk3_smoking_strata.py"
    spec = importlib.util.spec_from_file_location("dkk3_strata", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load DKK3 stratified effect script.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_dkk3_figure_module():
    path = ROOT / "scripts" / "09_figures" / "00_plot_gse292993_dkk3_compartments.py"
    spec = importlib.util.spec_from_file_location("dkk3_figure", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load DKK3 figure script.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_cell_composition_figure_module():
    path = ROOT / "scripts" / "09_figures" / "01_plot_gse292993_cell_composition_compartments.py"
    spec = importlib.util.spec_from_file_location("cell_composition_figure", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load cell-composition figure script.")
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
                payload = gzip.compress(dcc_text.encode("utf-8"))
                info = tarfile.TarInfo("nested/DSP-TEST-A01.dcc.gz")
                info.size = len(payload)
                archive.addfile(info, io.BytesIO(payload))

            pkc = geomx_module.inspect_pkc(pkc_path, "DKK3")
            dcc, rows = geomx_module.inspect_dcc_inputs(tar_path, root / "dcc", "DKK3")

        self.assertTrue(pkc["gene_present"])
        self.assertEqual(dcc["dcc_count"], 1)
        self.assertEqual(dcc["dccs_with_gene"], 1)
        self.assertEqual(dcc["raw_tar_extension_counts"][".dcc.gz"], 1)
        self.assertEqual(rows[0]["roi_id_guess"], "DSP-TEST-A01")
        self.assertTrue(rows[0]["contains_negative_probe_text"])

    def test_geomx_profile_parses_soft_and_merges_gsm(self):
        profile_module = load_geomx_profile_module()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            soft_path = root / "GSE292993_family.soft.gz"
            soft_text = "\n".join(
                [
                    "^SAMPLE = GSM8872219",
                    "!Sample_title = DSP-1001660011972-E-A01",
                    "!Sample_source_name_ch1 = lung",
                    "!Sample_characteristics_ch1 = diagnosis: COPD",
                    "!Sample_characteristics_ch1 = compartment: parenchyma",
                    "",
                ]
            )
            with gzip.open(soft_path, "wt", encoding="utf-8", newline="\n") as handle:
                handle.write(soft_text)

            dcc_path = root / "GSM8872219_DSP-1001660011972-E-A01.dcc.gz"
            with gzip.open(dcc_path, "wt", encoding="utf-8", newline="\n") as handle:
                handle.write("[FileHeader]\nSample_ID\tDSP-1001660011972-E-A01\n")

            geo_rows = profile_module.parse_soft_samples(soft_path)
            dcc_row = {
                "dcc_filename": dcc_path.name,
                "dcc_id": profile_module.dcc_id_from_name(dcc_path.name),
                "geo_accession": profile_module.gsm_from_name(dcc_path.name),
            }
            merged = profile_module.merge_dcc_and_geo([dcc_row], geo_rows)
            profile = profile_module.profile_dcc(dcc_path)

        self.assertEqual(geo_rows[0]["characteristics_diagnosis"], "COPD")
        self.assertTrue(merged[0]["metadata_matched"])
        self.assertEqual(merged[0]["characteristics_compartment"], "parenchyma")
        self.assertEqual(profile["key_values"]["Sample_ID"], "DSP-1001660011972-E-A01")

    def test_geomx_dcc_qc_resolves_pkc_codes_and_counts(self):
        dcc_qc_module = load_geomx_dcc_qc_module()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pkc_path = root / "panel.pkc.gz"
            pkc = {
                "Targets": [
                    {
                        "DisplayName": "DKK3",
                        "CodeClass": "Endogenous",
                        "RTS_ID": "RTS0020886",
                    },
                    {
                        "DisplayName": "ACTB",
                        "CodeClass": "Housekeeping",
                        "RTS_ID": "RTS0020894",
                    },
                ]
            }
            with gzip.open(pkc_path, "wt", encoding="utf-8", newline="\n") as handle:
                json.dump(pkc, handle)

            dcc_path = root / "GSM8872219_DSP-1001660011972-E-A01.dcc.gz"
            dcc_text = "\n".join(
                [
                    "<Scan_Attributes>",
                    "ID,DSP-1001660011972-E-A01",
                    "Plate_ID,1001660011972",
                    "Well,A01",
                    "</Scan_Attributes>",
                    "<NGS_Processing_Attributes>",
                    "Raw,10069",
                    "Trimmed,9905",
                    "Stitched,6315",
                    "Aligned,5995",
                    "umiQ30,0.9951",
                    "</NGS_Processing_Attributes>",
                    "<Code_Summary>",
                    "RTS0020886,17",
                    "RTS0020894,22",
                    "</Code_Summary>",
                    "",
                ]
            )
            with gzip.open(dcc_path, "wt", encoding="utf-8", newline="\n") as handle:
                handle.write(dcc_text)

            pkc_summary, pkc_rows = dcc_qc_module.parse_pkc_code_map(pkc_path, "DKK3")
            gene_codes = {
                row["code_id"] for row in pkc_rows if row["is_primary_gene"]
            }
            row, code_counts = dcc_qc_module.parse_dcc(dcc_path, gene_codes)

        self.assertEqual(pkc_summary["gene_code_count"], 1)
        self.assertEqual(gene_codes, {"RTS0020886"})
        self.assertEqual(code_counts["RTS0020886"], 17)
        self.assertEqual(row["primary_gene_counts"], 17)
        self.assertEqual(row["aligned_reads"], 5995)
        self.assertAlmostEqual(row["trimmed_fraction"], 9905 / 10069)
        self.assertAlmostEqual(row["stitched_fraction"], 6315 / 9905)
        self.assertAlmostEqual(row["aligned_fraction_stitched"], 5995 / 6315)

    def test_geomx_roi_qc_flags_low_quality_rows(self):
        roi_qc_module = load_geomx_roi_qc_module()

        class Args:
            min_aligned_reads = 100000
            min_code_counts = 10000
            min_total_code_counts = 10000
            min_trimmed_fraction = 0.90
            min_stitched_fraction = 0.80
            min_aligned_fraction_stitched = 0.80
            min_umi_q30 = 0.98
            min_rts_q30 = 0.98

        pass_row = {
            "metadata_matched": "True",
            "qc_metrics_matched": "True",
            "aligned_reads": "200000",
            "n_code_counts": "18000",
            "total_code_counts": "50000",
            "trimmed_fraction": "0.98",
            "stitched_fraction": "0.90",
            "aligned_fraction_stitched": "0.95",
            "umi_q30": "0.995",
            "rts_q30": "0.994",
        }
        fail_row = dict(pass_row)
        fail_row["aligned_reads"] = "1"
        fail_row["n_code_counts"] = "1"

        include, reasons = roi_qc_module.qc_flags(pass_row, Args)
        fail_include, fail_reasons = roi_qc_module.qc_flags(fail_row, Args)

        self.assertTrue(include)
        self.assertFalse(reasons)
        self.assertFalse(fail_include)
        self.assertIn("aligned_reads_below_min", fail_reasons)
        self.assertIn("n_code_counts_below_min", fail_reasons)
        self.assertEqual(roi_qc_module.disease_group("COPD"), "COPD")
        self.assertEqual(roi_qc_module.disease_group("Non Smoker"), "Control")
        self.assertEqual(roi_qc_module.disease_group("Smoker"), "Control")

    def test_geomx_loq_metrics_flags_background_detection(self):
        loq_module = load_geomx_loq_module()
        metrics = loq_module.loq_metrics(
            [1, 2, 3, 4],
            pseudocount=1.0,
            sd_multiplier=2.0,
        )
        rows = [
            {
                "include_qc": "True",
                "diagnosis_group": "COPD",
                "compartment_guess": "parenchyma",
                "dkk3_above_geometric_loq": True,
                "dkk3_above_arithmetic_loq": False,
            },
            {
                "include_qc": "False",
                "diagnosis_group": "COPD",
                "compartment_guess": "parenchyma",
                "dkk3_above_geometric_loq": True,
                "dkk3_above_arithmetic_loq": True,
            },
        ]
        summary = loq_module.summarize_by_group(
            rows, ["diagnosis_group", "compartment_guess"]
        )

        self.assertEqual(metrics["negative_probe_n"], 4)
        self.assertGreater(metrics["dkk3_geometric_loq"], 0)
        self.assertEqual(summary[0]["n_include_qc"], 1)
        self.assertEqual(summary[0]["n_dkk3_above_geometric_loq"], 1)
        self.assertEqual(summary[0]["n_dkk3_above_arithmetic_loq"], 0)

    def test_geomx_qc_plot_helpers_count_and_transform(self):
        plot_module = load_geomx_qc_plot_module()
        rows = [
            {
                "geo_accession": "GSM1",
                "diagnosis_guess": "COPD",
                "compartment_guess": "parenchyma",
                "donor_guess": "P1",
                "include_qc": "True",
                "total_code_counts": "99",
            },
            {
                "geo_accession": "GSM2",
                "diagnosis_guess": "Non Smoker",
                "compartment_guess": "unknown",
                "donor_guess": "P2",
                "include_qc": "False",
                "total_code_counts": "9",
            },
        ]
        metric = {
            "column": "total_code_counts",
            "transform": "log10p1",
        }
        grouped = plot_module.metric_values(rows, metric, "diagnosis_guess")
        counts = plot_module.count_table(
            rows, "diagnosis_guess", "compartment_guess"
        )
        donor_counts = plot_module.donor_count_table(
            rows, "diagnosis_guess", "compartment_guess"
        )
        summary = plot_module.qc_summary(rows)

        self.assertEqual(grouped["COPD"][0]["plot_value"], 2.0)
        self.assertEqual(counts["Non Smoker"]["unknown"], 1)
        self.assertEqual(donor_counts["COPD"]["parenchyma"], 1)
        self.assertEqual(
            plot_module.display_label("unknown", "diagnosis_guess"),
            "unknown diagnosis",
        )
        self.assertEqual(
            plot_module.display_label("unknown", "compartment_guess"),
            "unknown compartment",
        )
        self.assertEqual(summary["n_pass_qc"], 1)
        self.assertIsNone(
            summary["standard_ngs_like_metrics_available"]["mitochondrial_percent"]
        )

    def test_geomx_wta_matrix_feature_aggregation(self):
        matrix_module = load_geomx_matrix_module()
        pkc_rows = [
            {
                "code_id": "RTS1",
                "target": "DKK3",
                "code_class": "Endogenous",
            },
            {
                "code_id": "RTS2",
                "target": "DKK3",
                "code_class": "Endogenous",
            },
            {
                "code_id": "NEG1",
                "target": "NegProbe",
                "code_class": "Negative",
            },
        ]
        manifest, code_to_target, dropped = matrix_module.feature_manifest(
            pkc_rows, include_controls=False
        )
        counts = matrix_module.aggregate_counts(
            {"RTS1": 2, "RTS2": 3, "NEG1": 100}, code_to_target
        )
        logcpm = matrix_module.logcpm_matrix({"ROI1": counts}, ["DKK3"])

        self.assertEqual([row["target"] for row in manifest], ["DKK3"])
        self.assertEqual(dropped, [])
        self.assertEqual(counts["DKK3"], 5)
        self.assertGreater(logcpm["ROI1"]["DKK3"], 0)

    def test_geomx_wta_matrix_drops_implausibly_broad_pkc_targets(self):
        matrix_module = load_geomx_matrix_module()
        pkc_rows = [
            {
                "code_id": f"RTS{index}",
                "target": "Endogenous",
                "code_class": "Endogenous",
            }
            for index in range(30)
        ]
        manifest, code_to_target, dropped = matrix_module.feature_manifest(
            pkc_rows, include_controls=False, max_codes_per_target=20
        )

        self.assertEqual(manifest, [])
        self.assertEqual(code_to_target, {})
        self.assertEqual(dropped[0]["target"], "Endogenous")

    def test_scrna_reference_audit_helpers_resolve_overlap(self):
        audit_module = load_scrna_audit_module()

        class FakeSeries:
            def __init__(self, values):
                self.values = values

            def tolist(self):
                return self.values

        class FakeVar(dict):
            @property
            def columns(self):
                return list(self.keys())

        class FakeAdata:
            var_names = ["ENSG1", "ENSG2", "OTHER"]
            var = FakeVar({"gene_symbols": FakeSeries(["DKK3", "AGER", "BAD"])})

        self.assertEqual(
            audit_module.resolve_column(["Cell Type", "donor"], ["cell type"]),
            "Cell Type",
        )
        overlap = audit_module.select_gene_overlap(FakeAdata(), ["DKK3", "AGER"])

        self.assertEqual(overlap["source"], "var.gene_symbols")
        self.assertEqual(overlap["n_overlap"], 2)
        self.assertEqual(overlap["common_genes"], ["AGER", "DKK3"])

    def test_scrna_reference_audit_reads_large_geomx_manifest_fields(self):
        audit_module = load_scrna_audit_module()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "feature_manifest.csv"
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "target",
                        "n_codes",
                        "code_ids",
                        "is_control_feature",
                    ],
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "target": "DKK3",
                        "n_codes": "1",
                        "code_ids": "RTS1",
                        "is_control_feature": "False",
                    }
                )
                writer.writerow(
                    {
                        "target": "BROAD_RECORD",
                        "n_codes": "999",
                        "code_ids": ";".join(f"RTS{index}" for index in range(50000)),
                        "is_control_feature": "False",
                    }
                )

            genes = audit_module.read_geomx_genes(path)

        self.assertEqual(genes, ["DKK3"])

    def test_parenchyma_deconvolution_helpers_filter_and_correlate(self):
        deconv_module = load_deconv_module()
        rows = [
            {
                "geo_accession": "GSM1",
                "include_qc": "True",
                "diagnosis_group": "COPD",
                "diagnosis_guess": "COPD",
                "compartment_guess": "parenchyma",
                "donor_guess": "P1",
                "fraction_fibroblast": "0.8",
                "log1p_dkk3_cpm": "4.0",
            },
            {
                "geo_accession": "GSM2",
                "include_qc": "True",
                "diagnosis_group": "Control",
                "diagnosis_guess": "Non Smoker",
                "compartment_guess": "airway",
                "donor_guess": "P2",
                "fraction_fibroblast": "0.1",
                "log1p_dkk3_cpm": "1.0",
            },
            {
                "geo_accession": "GSM3",
                "include_qc": "False",
                "diagnosis_group": "COPD",
                "diagnosis_guess": "COPD",
                "compartment_guess": "parenchyma",
                "donor_guess": "P3",
                "fraction_fibroblast": "0.2",
                "log1p_dkk3_cpm": "2.0",
            },
        ]
        metadata = deconv_module.roi_metadata(rows, "parenchyma")
        donor_rows = deconv_module.donor_summaries(
            [rows[0]], ["fraction_fibroblast"]
        )
        correlations = deconv_module.correlation_rows(
            [
                {"fraction_fibroblast": "0.1", "log1p_dkk3_cpm": "1.0"},
                {"fraction_fibroblast": "0.5", "log1p_dkk3_cpm": "2.0"},
                {"fraction_fibroblast": "0.9", "log1p_dkk3_cpm": "3.0"},
            ],
            ["fraction_fibroblast"],
            unit="roi",
        )

        self.assertEqual(deconv_module.safe_name("AT2 / epithelial"), "at2_epithelial")
        self.assertEqual(list(metadata), ["GSM1"])
        self.assertEqual(donor_rows[0]["median_fraction_fibroblast"], 0.8)
        self.assertAlmostEqual(correlations[0]["spearman_rho"], 1.0)

    def test_scanpy_notebook_inspection_ignores_outputs_and_finds_annotation_clues(self):
        notebook_module = load_scanpy_notebook_inspection_module()
        notebook = {
            "cells": [
                {
                    "cell_type": "code",
                    "source": [
                        "adata = sc.read_10x_h5('GSM1_filtered_feature_bc_matrix.h5')\n",
                    ],
                    "outputs": [
                        {
                            "data": {
                                "image/png": "cell_type marker_gene leiden read_h5ad"
                            }
                        }
                    ],
                },
                {
                    "cell_type": "code",
                    "source": [
                        "marker_gene_dict = {'fibroblast': ['COL1A1']}\n",
                        "adata.obs['cell_type'] = adata.obs['leiden'].map(labels)\n",
                    ],
                },
                {
                    "cell_type": "code",
                    "source": [
                        "adata.write_h5ad('atlas_full_annotated.h5ad')\n",
                        "with open('output/adata_harmony_annotated_cr8', 'wb') as file:\n",
                        "    pickle.dump(adata, file)\n",
                    ],
                },
            ]
        }
        with tempfile.TemporaryDirectory() as directory:
            zip_path = Path(directory) / "scanpy_workflow.zip"
            with zipfile.ZipFile(zip_path, mode="w") as archive:
                archive.writestr("2_celltype_annotation.ipynb", json.dumps(notebook))
            overall, summaries, matches, artifacts = notebook_module.summarize_archive(
                zip_path, notebook_module.DEFAULT_TERMS
            )

        self.assertEqual(overall["n_notebooks"], 1)
        self.assertEqual(
            overall["notebooks_with_write_h5ad_or_write"],
            ["2_celltype_annotation.ipynb"],
        )
        self.assertTrue(summaries[0]["has_celltype_or_annotation_terms"])
        self.assertIn("atlas_full_annotated.h5ad", summaries[0]["mentioned_paths"])
        self.assertIn("output/adata_harmony_annotated_cr8", summaries[0]["written_artifact_paths"])
        self.assertEqual(artifacts[0]["direction"], "write")
        self.assertTrue(any("marker_gene_dict" in row["relevant_lines"] for row in matches))

    def test_annotation_code_export_writes_selected_notebook_code(self):
        export_module = load_annotation_code_export_module()
        notebook = {
            "cells": [
                {
                    "cell_type": "code",
                    "source": [
                        "import scvi\n",
                        "scvi.settings.seed = 0\n",
                        "import mudata as mu\n",
                        "!rm -r \"input/data_cellranger8/.DS_Store\"\n",
                        "celldict_level1 = {'Parenchyma': ['0', '1']}\n",
                        "adata.obs['celltype_level1'] = 'Parenchyma'\n",
                        "with open('output/adata_harmony_annotated_cr8', 'wb') as file:\n",
                        "    pickle.dump(adata, file)\n",
                    ],
                }
            ]
        }
        with tempfile.TemporaryDirectory() as directory:
            directory_path = Path(directory)
            zip_path = directory_path / "scanpy_workflow.zip"
            output_dir = directory_path / "code"
            with zipfile.ZipFile(zip_path, mode="w") as archive:
                archive.writestr("scanpy_workflow/2_celltype_annotation.ipynb", json.dumps(notebook))
                archive.writestr("scanpy_workflow/unrelated.ipynb", json.dumps(notebook))
            summary, rows = export_module.export_annotation_code(
                zip_path,
                output_dir,
                ["2_celltype_annotation"],
            )

            code_file = Path(summary["code_files"][0])
            exported = code_file.read_text(encoding="utf-8")
            py_compile.compile(str(code_file), doraise=True)

        self.assertEqual(summary["n_selected_notebooks"], 1)
        self.assertEqual(summary["cells_with_celldict"], 1)
        self.assertEqual(summary["cells_with_open_output"], 1)
        self.assertEqual(summary["n_notebook_only_lines_commented"], 4)
        self.assertTrue(rows[0]["has_celltype_assignment"])
        self.assertIn("# %% [cell 0]", exported)
        self.assertIn("pass  # NOTEBOOK-OPTIONAL-SCVI", exported)
        self.assertIn("pass  # NOTEBOOK-OPTIONAL-MUDATA", exported)
        self.assertIn("pass  # NOTEBOOK-ONLY: !rm", exported)
        self.assertIn("celldict_level1", exported)

    def test_author_celltype_map_parser_tracks_reused_celldict_contexts(self):
        map_module = load_author_celltype_map_module()
        text = """
celldict_level1 = {
    'Parenchyma': ['0', '1'],
    'Endothelial': ['2'],
}
adata.obs['celltype_level1']=np.NaN
for i in celldict_level1.keys():
    pass

celldict_level1 = {
    'AT1': ['0'],
    'AT2': ['1'],
    'Fibroblast': ['2', '3'],
}
parenchyma.obs['parenchyma_celltype_level1']=np.NaN
for i in celldict_level1.keys():
    pass
"""
        rows, blocks = map_module.extract_celltype_maps_from_text(
            text,
            source="synthetic.py",
        )
        parenchyma_rows = [
            row
            for row in rows
            if row["assigned_obs_column"] == "parenchyma_celltype_level1"
        ]

        self.assertEqual(len(blocks), 2)
        self.assertEqual(blocks[0]["assigned_object"], "adata")
        self.assertEqual(blocks[1]["assigned_object"], "parenchyma")
        self.assertEqual(len(parenchyma_rows), 3)
        self.assertEqual(
            next(row for row in parenchyma_rows if row["celltype_label"] == "Fibroblast")[
                "cluster_ids"
            ],
            "2;3",
        )

    def test_author_signature_builder_validates_labels_and_anndata_like_objects(self):
        signature_module = load_author_signature_module()

        class FakeAdata:
            obs = {"celltype": ["AT1"]}
            var = {"gene": ["AGER"]}
            X = [[1]]
            n_obs = 1

        fake = FakeAdata()

        self.assertTrue(signature_module.valid_label("AT1"))
        self.assertFalse(signature_module.valid_label("nan"))
        self.assertFalse(signature_module.valid_label(""))
        self.assertTrue(signature_module.has_anndata_interface(fake))
        self.assertIs(signature_module.coerce_anndata_like({"adata": fake}), fake)
        with self.assertRaises(TypeError):
            signature_module.coerce_anndata_like({"not_adata": object()})

    def test_dkk3_summary_enriches_and_counts_donors(self):
        dkk3_module = load_dkk3_summary_module()
        rows = dkk3_module.enrich_roi_rows(
            [
                {
                    "include_qc": "True",
                    "diagnosis_group": "COPD",
                    "compartment_guess": "parenchyma",
                    "donor_guess": "P1",
                    "dkk3_count": "10",
                    "total_code_counts": "1000",
                    "dkk3_above_geometric_loq": "True",
                    "dkk3_above_arithmetic_loq": "False",
                },
                {
                    "include_qc": "False",
                    "diagnosis_group": "COPD",
                    "compartment_guess": "parenchyma",
                    "donor_guess": "P1",
                    "dkk3_count": "100",
                    "total_code_counts": "1000",
                    "dkk3_above_geometric_loq": "True",
                    "dkk3_above_arithmetic_loq": "True",
                },
                {
                    "include_qc": "True",
                    "diagnosis_group": "Control",
                    "compartment_guess": "airway",
                    "donor_guess": "P2",
                    "dkk3_count": "0",
                    "total_code_counts": "1000",
                    "dkk3_above_geometric_loq": "False",
                    "dkk3_above_arithmetic_loq": "False",
                },
                {
                    "include_qc": "True",
                    "diagnosis_group": "Control",
                    "compartment_guess": "unknown",
                    "donor_guess": "P3",
                    "dkk3_count": "5",
                    "total_code_counts": "1000",
                    "dkk3_above_geometric_loq": "False",
                    "dkk3_above_arithmetic_loq": "False",
                },
            ]
        )
        primary = dkk3_module.primary_rows(rows)
        summary = dkk3_module.summarize_group(
            primary, ["diagnosis_group", "compartment_guess"]
        )
        donor_counts = dkk3_module.donor_counts_by_column(primary, "diagnosis_group")

        self.assertEqual(len(primary), 2)
        self.assertAlmostEqual(rows[0]["dkk3_cpm"], 10000)
        copd_summary = next(row for row in summary if row["diagnosis_group"] == "COPD")
        self.assertEqual(copd_summary["n_dkk3_above_geometric_loq"], 1)
        self.assertEqual(
            donor_counts,
            [
                {"diagnosis_group": "COPD", "n_donors": 1},
                {"diagnosis_group": "Control", "n_donors": 1},
            ],
        )

    def test_dkk3_effects_compare_donor_groups(self):
        effects_module = load_dkk3_effect_module()
        rng = random.Random(123)
        rows = [
            {
                "donor_guess": "C1",
                "diagnosis_group": "COPD",
                "compartment_guess": "parenchyma",
                "median_log1p_dkk3_cpm": "3.0",
            },
            {
                "donor_guess": "C2",
                "diagnosis_group": "COPD",
                "compartment_guess": "parenchyma",
                "median_log1p_dkk3_cpm": "4.0",
            },
            {
                "donor_guess": "N1",
                "diagnosis_group": "Control",
                "compartment_guess": "parenchyma",
                "median_log1p_dkk3_cpm": "1.0",
            },
            {
                "donor_guess": "N2",
                "diagnosis_group": "Control",
                "compartment_guess": "parenchyma",
                "median_log1p_dkk3_cpm": "2.0",
            },
        ]
        result = effects_module.test_compartment_metric(
            rows,
            "parenchyma",
            "median_log1p_dkk3_cpm",
            permutations=100,
            bootstraps=100,
            rng=rng,
        )

        self.assertEqual(result["n_copd_donors"], 2)
        self.assertEqual(result["n_control_donors"], 2)
        self.assertEqual(result["mean_difference_copd_minus_control"], 2.0)
        self.assertIsNotNone(result["permutation_p_two_sided"])

    def test_dkk3_strata_pairwise_comparison(self):
        strata_module = load_dkk3_strata_module()
        effects_module = load_dkk3_effect_module()
        rng = random.Random(456)
        rows = [
            {
                "donor_guess": "C1",
                "diagnosis_guess": "COPD",
                "compartment_guess": "parenchyma",
                "median_log1p_dkk3_cpm": "4.0",
            },
            {
                "donor_guess": "C2",
                "diagnosis_guess": "COPD",
                "compartment_guess": "parenchyma",
                "median_log1p_dkk3_cpm": "5.0",
            },
            {
                "donor_guess": "S1",
                "diagnosis_guess": "Smoker",
                "compartment_guess": "parenchyma",
                "median_log1p_dkk3_cpm": "3.0",
            },
            {
                "donor_guess": "N1",
                "diagnosis_guess": "Non Smoker",
                "compartment_guess": "parenchyma",
                "median_log1p_dkk3_cpm": "1.0",
            },
            {
                "donor_guess": "N2",
                "diagnosis_guess": "Non Smoker",
                "compartment_guess": "parenchyma",
                "median_log1p_dkk3_cpm": "2.0",
            },
        ]
        result = strata_module.test_pair_metric(
            rows,
            label_column="diagnosis_guess",
            label_a="COPD",
            label_b="Non Smoker",
            compartment="parenchyma",
            metric="median_log1p_dkk3_cpm",
            permutations=100,
            bootstraps=100,
            rng=rng,
            effects=effects_module,
        )

        self.assertEqual(strata_module.parse_pairs(None, ["COPD", "Non Smoker", "Smoker"])[0], ("COPD", "Non Smoker"))
        self.assertEqual(result["n_label_a_donors"], 2)
        self.assertEqual(result["n_label_b_donors"], 2)
        self.assertEqual(result["mean_difference_label_a_minus_label_b"], 3.0)

    def test_dkk3_figure_helpers_filter_and_format(self):
        figure_module = load_dkk3_figure_module()
        donor_rows = [
            {
                "donor_guess": "P1",
                "diagnosis_guess": "COPD",
                "compartment_guess": "parenchyma",
                "median_log1p_dkk3_cpm": "4.0",
            },
            {
                "donor_guess": "P2",
                "diagnosis_guess": "Smoker",
                "compartment_guess": "airway",
                "median_log1p_dkk3_cpm": "5.0",
            },
            {
                "donor_guess": "P3",
                "diagnosis_guess": "Non Smoker",
                "compartment_guess": "parenchyma",
                "median_log1p_dkk3_cpm": "3.0",
            },
            {
                "donor_guess": "P4",
                "diagnosis_guess": "COPD",
                "compartment_guess": "vessel",
                "median_log1p_dkk3_cpm": "2.0",
            },
        ]
        effect_rows = [
            {
                "compartment_guess": "parenchyma",
                "metric": "median_log1p_dkk3_cpm",
                "label_a": "COPD",
                "label_b": "Non Smoker",
                "mean_difference_label_a_minus_label_b": "1.0",
                "permutation_p_two_sided": "0.031",
            }
        ]
        grouped = figure_module.compartment_values(
            donor_rows, "median_log1p_dkk3_cpm", "parenchyma"
        )
        effects = figure_module.effect_lookup(
            effect_rows, "parenchyma", "median_log1p_dkk3_cpm"
        )

        self.assertEqual(len(grouped["COPD"]), 1)
        self.assertEqual(len(grouped["Smoker"]), 0)
        self.assertEqual(
            figure_module.parse_compartments("all"),
            ["airway", "parenchyma", "vessel"],
        )
        self.assertEqual(
            figure_module.parse_compartments("diagnostic"),
            ["airway", "parenchyma", "vessel", "unknown"],
        )
        self.assertEqual(
            figure_module.parse_compartments("parenchyma,vessel"),
            ["parenchyma", "vessel"],
        )
        self.assertEqual(
            figure_module.output_stem(["airway", "parenchyma", "vessel"]),
            "gse292993_dkk3_all_compartments_donor_signal",
        )
        self.assertEqual(
            figure_module.output_stem(["airway", "parenchyma", "vessel", "unknown"]),
            "gse292993_dkk3_all_compartments_with_unknown_donor_signal",
        )
        plot_summary = figure_module.plot_count_summary(
            donor_rows, ["airway", "parenchyma", "vessel"]
        )
        self.assertEqual(plot_summary["missing_requested_compartments"], [])
        self.assertIn(
            "vessel", plot_summary["available_compartments_in_donor_table"]
        )
        self.assertIn("COPD-Non Smoker", figure_module.metric_effect_text(effects))
        self.assertEqual(
            figure_module.stable_jitter("P1"),
            figure_module.stable_jitter("P1"),
        )

    def test_cell_composition_helpers_make_donor_balanced_summary(self):
        figure_module = load_cell_composition_figure_module()
        manifest = [
            {"cell_type": "AT1", "fraction_column": "fraction_at1"},
            {"cell_type": "Fibroblast", "fraction_column": "fraction_fibroblast"},
        ]
        roi_rows = [
            {
                "donor_guess": "D1",
                "diagnosis_guess": "COPD",
                "fraction_at1": "0.2",
                "fraction_fibroblast": "0.8",
            },
            {
                "donor_guess": "D1",
                "diagnosis_guess": "COPD",
                "fraction_at1": "0.4",
                "fraction_fibroblast": "0.6",
            },
            {
                "donor_guess": "D2",
                "diagnosis_guess": "Non Smoker",
                "fraction_at1": "0.9",
                "fraction_fibroblast": "0.1",
            },
            {
                "donor_guess": "D3",
                "diagnosis_guess": "unknown",
                "fraction_at1": "0.5",
                "fraction_fibroblast": "0.5",
            },
        ]
        donor_rows = figure_module.donor_fraction_rows(
            roi_rows=roi_rows,
            manifest_rows=manifest,
            compartment="parenchyma",
        )
        summary = figure_module.composition_summary_rows(
            donor_rows,
            manifest_by_compartment={"parenchyma": manifest},
            compartments=["parenchyma"],
            labels=["Non Smoker", "COPD"],
        )
        lookup = figure_module.value_lookup(summary)
        stacked = figure_module.stack_values(
            lookup,
            compartment="parenchyma",
            diagnosis="COPD",
            cell_types=["AT1", "Fibroblast"],
            normalize=True,
        )

        self.assertEqual(figure_module.parse_compartments("all"), ["airway", "parenchyma", "vessel"])
        self.assertEqual(len(donor_rows), 2)
        self.assertAlmostEqual(
            next(row for row in donor_rows if row["donor_guess"] == "D1")[
                "mean_fraction_at1"
            ],
            0.3,
        )
        self.assertAlmostEqual(
            lookup[("parenchyma", "COPD", "Fibroblast")],
            0.7,
        )
        self.assertAlmostEqual(sum(stacked), 1.0)
        self.assertEqual(
            figure_module.celltype_order({"parenchyma": manifest}),
            ["AT1", "Fibroblast"],
        )

    def test_author_input_helper_identifies_ribosomal_prefixes(self):
        helper_module = load_author_input_helper_module()
        pattern = helper_module.ribosomal_gene_pattern(("RPL", "RPS", "MRPL", "MRPS"))

        self.assertTrue(pattern.match("RPL13A"))
        self.assertTrue(pattern.match("RPS27"))
        self.assertTrue(pattern.match("MRPL12"))
        self.assertTrue(pattern.match("MRPS18B"))
        self.assertFalse(pattern.match("DKK3"))
        self.assertFalse(pattern.match("RPRD1A"))

    def test_author_input_helper_finds_zip_sidecar_by_basename(self):
        helper_module = load_author_input_helper_module()
        with tempfile.TemporaryDirectory() as tmpdir:
            zip_path = Path(tmpdir) / "scanpy_workflow.zip"
            with zipfile.ZipFile(zip_path, "w") as archive:
                archive.writestr(
                    "scanpy_workflow/input/meta_cr8.csv",
                    "sample,group\ns1,COPD\n",
                )

            self.assertEqual(
                helper_module.find_zip_member(zip_path, "meta_cr8.csv"),
                "scanpy_workflow/input/meta_cr8.csv",
            )
            self.assertIsNone(helper_module.find_zip_member(zip_path, "missing.csv"))

    def test_author_input_helper_discovers_notebook_sidecar_reads(self):
        helper_module = load_author_input_helper_module()
        notebook = {
            "cells": [
                {
                    "cell_type": "code",
                    "source": [
                        "meta = pd.read_csv('input/meta_cr8.csv')\n",
                        "ribo = pd.read_table(\"input/GOCC_RIBOSOMAL_SUBUNIT.v2023.1.Hs.csv\", header=None)\n",
                        "with open('output/adata_harmony_integrated_cr8', 'wb') as handle:\n",
                        "    pass\n",
                    ],
                }
            ]
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            zip_path = Path(tmpdir) / "scanpy_workflow.zip"
            with zipfile.ZipFile(zip_path, "w") as archive:
                archive.writestr("scanpy_workflow/1_preprocessing.ipynb", json.dumps(notebook))
                archive.writestr("scanpy_workflow/input/meta_cr8.csv", "sample,group\ns1,COPD\n")

            refs = helper_module.discover_notebook_file_references(zip_path)
            read_sidecars = helper_module.unique_read_sidecars(refs)

            self.assertIn("input/meta_cr8.csv", read_sidecars)
            self.assertIn(
                "input/GOCC_RIBOSOMAL_SUBUNIT.v2023.1.Hs.csv",
                read_sidecars,
            )
            self.assertNotIn("output/adata_harmony_integrated_cr8", read_sidecars)
            self.assertEqual(
                helper_module.find_zip_member_for_path(zip_path, "input/meta_cr8.csv"),
                "scanpy_workflow/input/meta_cr8.csv",
            )


if __name__ == "__main__":
    unittest.main()
