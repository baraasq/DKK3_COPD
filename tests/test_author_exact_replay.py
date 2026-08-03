from __future__ import annotations

import csv
import gzip
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


def load_script(name: str, filename: str):
    path = ROOT / "scripts" / "06_celltype" / filename
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


PREPARE = load_script(
    "gse302339_author_exact_prepare_test",
    "15_prepare_gse302339_author_exact_replay.py",
)
EXPORT = load_script(
    "gse302339_author_exact_export_test",
    "16_export_gse302339_author_exact_annotations.py",
)


def code_cell(source: str = "", output: str = "") -> dict:
    outputs = [{"output_type": "stream", "text": [output]}] if output else []
    return {
        "cell_type": "code",
        "source": [source],
        "outputs": outputs,
    }


def synthetic_notebooks(n_samples: int = 65) -> tuple[dict, dict]:
    post_doublet = [1000 + index for index in range(n_samples - 1)]
    post_doublet.append(160620 - sum(post_doublet))
    raw = [value + 100 for value in post_doublet]
    post_low_quality = [value + 50 for value in post_doublet]

    sample_output = "\n".join(
        (
            f"Total number of cells: {raw_count}\n"
            f"Number of cells after filtering of low quality cells: {low_quality}\n"
            f"({final_count}, {13000 + position})"
        )
        for position, (raw_count, low_quality, final_count) in enumerate(
            zip(raw, post_low_quality, post_doublet, strict=True)
        )
    )
    preprocessing = {
        "metadata": {"language_info": {"version": "3.12.8"}},
        "cells": [
            code_cell(),
            code_cell(
                output=(
                    "scanpy==1.10.4 anndata==0.11.1 numpy==1.26.4 "
                    "pandas==2.2.3"
                )
            ),
            code_cell(),
            code_cell(output=sample_output),
            code_cell(output="(160620, 18941)\n(160620, 2323)"),
        ],
    }
    annotation = {
        "cells": [
            code_cell(
                output=(
                    "finished: found 41 clusters and added\n"
                    "finished: found 38 clusters and added"
                )
            )
        ]
    }
    return preprocessing, annotation


class FakeSeries:
    def __init__(self, values: list[str]):
        self.values = values

    def astype(self, _dtype):
        return self

    def nunique(self) -> int:
        return len(set(self.values))

    def __iter__(self):
        return iter(self.values)


class FakeObs:
    def __init__(self, rows: list[tuple[str, dict]], series: dict[str, list[str]]):
        self._rows = rows
        self._series = series
        row_columns = {key for _, row in rows for key in row}
        self.columns = row_columns | set(series)

    def __getitem__(self, column: str) -> FakeSeries:
        return FakeSeries(self._series[column])

    def iterrows(self):
        return iter(self._rows)


class FakeAdata:
    def __init__(
        self,
        rows: list[tuple[str, dict]],
        cluster_ids: list[str],
        shape: tuple[int, int],
    ):
        self.obs = FakeObs(rows, {"leiden": cluster_ids})
        self.shape = shape
        self.n_obs = len(rows)


class NotebookContractTests(unittest.TestCase):
    def test_deposited_archive_contract_when_archive_is_available(self):
        zip_path = PREPARE.resolve_zip(None)
        if not zip_path.exists():
            self.skipTest("deposited Scanpy notebook archive is not present")

        preprocessing = PREPARE.read_notebook(
            zip_path, PREPARE.PREPROCESSING_MEMBER
        )
        annotation = PREPARE.read_notebook(zip_path, PREPARE.ANNOTATION_MEMBER)
        contract = PREPARE.notebook_contract(preprocessing, annotation)

        self.assertEqual(contract["n_author_samples"], 65)
        self.assertEqual(
            contract["raw_cell_totals_in_author_order"][:2], [2437, 3603]
        )
        self.assertEqual(contract["raw_cell_totals_in_author_order"][-1], 1450)
        self.assertEqual(
            contract["post_low_quality_cell_counts_in_author_order"][:2],
            [2043, 3098],
        )
        self.assertEqual(
            contract["post_doublet_shapes_in_author_order"][:2],
            [[1363, 13394], [2640, 17109]],
        )
        self.assertEqual(
            contract["post_doublet_shapes_in_author_order"][-1], [982, 14036]
        )
        self.assertEqual(
            sum(row[0] for row in contract["post_doublet_shapes_in_author_order"]),
            160620,
        )
        self.assertEqual(contract["full_cluster_ids"], [str(i) for i in range(62)])
        self.assertEqual(contract["parenchyma_cluster_count"], 41)
        self.assertEqual(contract["immune_cluster_count"], 38)

    def test_notebook_contract_recovers_all_author_checkpoints(self):
        preprocessing, annotation = synthetic_notebooks()
        contract = PREPARE.notebook_contract(preprocessing, annotation)

        self.assertEqual(contract["python"], "3.12.8")
        self.assertEqual(contract["n_author_samples"], 65)
        self.assertEqual(len(contract["raw_cell_totals_in_author_order"]), 65)
        self.assertEqual(
            len(contract["post_low_quality_cell_counts_in_author_order"]), 65
        )
        self.assertEqual(len(contract["post_doublet_shapes_in_author_order"]), 65)
        self.assertTrue(contract["raw_cell_totals_are_unique"])
        self.assertEqual(
            sum(row[0] for row in contract["post_doublet_shapes_in_author_order"]),
            160620,
        )
        self.assertEqual(contract["concat_shape"], [160620, 18941])
        self.assertEqual(contract["post_hvg_shape"], [160620, 2323])
        self.assertEqual(contract["full_cluster_count"], 62)
        self.assertEqual(contract["full_cluster_ids"], [str(i) for i in range(62)])
        self.assertEqual(contract["parenchyma_cluster_count"], 41)
        self.assertEqual(contract["immune_cluster_count"], 38)
        self.assertEqual(contract["harmony_iterations_in_notebook_log"], 26)

    def test_notebook_contract_rejects_incomplete_sample_output(self):
        preprocessing, annotation = synthetic_notebooks(n_samples=64)
        with self.assertRaisesRegex(ValueError, "Expected 65"):
            PREPARE.notebook_contract(preprocessing, annotation)

    def test_extracts_exact_185_gene_ribosomal_tuple(self):
        genes = [f"RIBO_{index:03d}" for index in range(185)]
        code = (
            "unrelated = ('NOT_A_RIBOSOMAL_LIST',)\n"
            "ribosome = adata.var_names.str.startswith("
            + repr(tuple(genes))
            + ")\n"
        )

        observed = PREPARE.extract_ribosomal_genes(code)

        self.assertEqual(observed, genes)
        self.assertEqual(len(observed), 185)


class ExportContractTests(unittest.TestCase):
    def test_exact_cluster_constants_and_sample_barcode_key(self):
        self.assertEqual(EXPORT.EXPECTED_SHAPE, (160620, 2323))
        self.assertEqual(EXPORT.EXPECTED_FULL_CLUSTERS, 62)
        self.assertEqual(EXPORT.EXPECTED_PARENCHYMA_CLUSTERS, 41)
        self.assertEqual(EXPORT.EXPECTED_IMMUNE_CLUSTERS, 38)
        self.assertEqual(EXPORT.raw_barcode("AACCGGTT-1-4"), "AACCGGTT-1")
        self.assertEqual(
            EXPORT.sample_barcode_key(
                "AACCGGTT-1-4",
                {"batch": "input/data_cellranger8/GSM1_sample"},
            ),
            "input/data_cellranger8/GSM1_sample::AACCGGTT-1",
        )

    def test_export_allows_author_megakarycyte_as_broad_only_label(self):
        parenchyma_labels = sorted(EXPORT.REQUIRED_PARENCHYMA_LABELS)
        full_rows: list[tuple[str, dict]] = []
        parenchyma_rows: list[tuple[str, dict]] = []
        for index, label in enumerate(parenchyma_labels):
            batch = f"input/data_cellranger8/parenchyma_{index}"
            cell_id = "AAAA-1"
            full_rows.append(
                (
                    cell_id,
                    {"batch": batch, "celltype_level1": "Parenchyma"},
                )
            )
            parenchyma_rows.append(
                (
                    cell_id,
                    {"batch": batch, "parenchyma_celltype_level1": label},
                )
            )

        immune_batch = "input/data_cellranger8/immune_0"
        megakarycyte_batch = "input/data_cellranger8/megakarycyte_0"
        full_rows.extend(
            [
                (
                    "CCCC-1",
                    {"batch": immune_batch, "celltype_level1": "Immune"},
                ),
                (
                    "GGGG-1",
                    {
                        "batch": megakarycyte_batch,
                        "celltype_level1": "Megakarycyte",
                    },
                ),
            ]
        )
        immune_rows = [
            (
                "CCCC-1",
                {
                    "batch": immune_batch,
                    "immune_celltype_level1": "Macrophage",
                },
            )
        ]

        objects = {
            "adata_harmony_annotated_cr8": FakeAdata(
                full_rows,
                [str(i) for i in range(62)],
                EXPORT.EXPECTED_SHAPE,
            ),
            "parenchyma_harmony_annotated_cr8": FakeAdata(
                parenchyma_rows,
                [str(i) for i in range(41)],
                (len(parenchyma_rows), EXPORT.EXPECTED_SHAPE[1]),
            ),
            "immune_harmony_annotated_cr8": FakeAdata(
                immune_rows,
                [str(i) for i in range(38)],
                (len(immune_rows), EXPORT.EXPECTED_SHAPE[1]),
            ),
        }

        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary = Path(temporary_directory)
            output_dir = temporary / "objects"
            output_dir.mkdir()
            for filename in objects:
                (output_dir / filename).touch()
            annotation_output = temporary / "annotations.tsv.gz"
            preprocessing_log = temporary / "preprocessing.log"
            preprocessing_log.write_text(
                "harmonypy - INFO - Converged after 26 iterations\n",
                encoding="utf-8",
            )
            annotation_log = temporary / "annotation.log"
            annotation_log.write_text(
                "Exact author major-lineage annotations completed and checkpointed.\n",
                encoding="utf-8",
            )
            results_meta = temporary / "results" / "meta"
            results_tables = temporary / "results" / "tables"
            results_meta.mkdir(parents=True)
            results_tables.mkdir(parents=True)

            def fake_load_pickle(path: Path):
                return objects[path.name]

            argv = [
                "16_export_gse302339_author_exact_annotations.py",
                "--output-dir",
                str(output_dir),
                "--annotation-output",
                str(annotation_output),
                "--preprocessing-log",
                str(preprocessing_log),
                "--annotation-log",
                str(annotation_log),
                "--strict",
            ]
            with (
                mock.patch.object(sys, "argv", argv),
                mock.patch.object(EXPORT, "project_path", side_effect=Path),
                mock.patch.object(EXPORT, "load_pickle", side_effect=fake_load_pickle),
                mock.patch.object(EXPORT, "load_config", return_value={}),
                mock.patch.object(
                    EXPORT,
                    "ensure_results_dirs",
                    return_value={"meta": results_meta, "tables": results_tables},
                ),
            ):
                self.assertEqual(EXPORT.main(), 0)

            with gzip.open(annotation_output, "rt", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            megakarycyte = next(
                row for row in rows if row["celltype_level1"] == "Megakarycyte"
            )
            self.assertEqual(megakarycyte["parenchyma_celltype_level1"], "")
            self.assertEqual(megakarycyte["immune_celltype_level1"], "")
            self.assertEqual(megakarycyte["author_celltype"], "Megakarycyte")


if __name__ == "__main__":
    unittest.main()
