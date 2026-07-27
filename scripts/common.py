from __future__ import annotations

import os
import tomllib
from pathlib import Path
from typing import Iterable, Mapping, Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = PROJECT_ROOT / "config" / "project.toml"


def load_config(path: Path = CONFIG_PATH) -> dict:
    with path.open("rb") as handle:
        return tomllib.load(handle)


def project_path(value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else PROJECT_ROOT / path


def configured_path(
    config: Mapping,
    key: str,
    *,
    environment_variable: str | None = None,
) -> Path:
    if environment_variable and os.environ.get(environment_variable):
        return Path(os.environ[environment_variable]).expanduser()
    return project_path(config["paths"][key])


def ensure_results_dirs(config: Mapping) -> dict[str, Path]:
    root = configured_path(config, "results_dir")
    paths = {
        "root": root,
        "meta": root / "meta",
        "tables": root / "tables",
        "figures": root / "figures",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def resolve_column(
    columns: Iterable[str],
    candidates: Sequence[str],
    *,
    label: str,
    required: bool = True,
) -> str | None:
    available = list(columns)
    exact = {str(column): str(column) for column in available}
    folded = {str(column).casefold(): str(column) for column in available}

    for candidate in candidates:
        if candidate in exact:
            return exact[candidate]
        if candidate.casefold() in folded:
            return folded[candidate.casefold()]

    if required:
        raise KeyError(
            f"Could not resolve {label}. Tried {list(candidates)}. "
            f"Available columns: {available}"
        )
    return None


def normalize_condition(value: object, config: Mapping) -> str:
    text = str(value).strip()
    if text.casefold() == str(config["project"]["disease_label"]).casefold():
        return config["project"]["disease_label"]
    control_labels = config["project"]["control_labels"]
    if text.casefold() in {str(label).casefold() for label in control_labels}:
        return "Control"
    return text


def expression_vector(matrix) -> "object":
    """Return a one-dimensional NumPy array from dense or sparse input."""
    import numpy as np

    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    return np.asarray(matrix).reshape(-1)


def scrinshot_cell_map_members(names: Iterable[str]) -> list[str]:
    return sorted(
        name
        for name in names
        if "Cell maps and histology/" in name
        and name.casefold().endswith(".csv")
    )
