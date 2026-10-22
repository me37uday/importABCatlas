from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def _description_version() -> str:
    for line in (ROOT / "DESCRIPTION").read_text(encoding="utf-8").splitlines():
        if line.startswith("Version:"):
            return line.split(":", 1)[1].strip()
    raise AssertionError("DESCRIPTION has no Version field")


def _backend_version() -> str:
    path = ROOT / "inst" / "python" / "abc_backend.py"
    spec = importlib.util.spec_from_file_location("iab_backend_for_version", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return str(module.VERSION)


def test_description_and_backend_versions_match() -> None:
    assert _description_version() == "0.2.0.9007"
    assert _backend_version() == _description_version()


def test_source_tree_contains_no_installed_lazyload_databases() -> None:
    forbidden = [
        path
        for path in ROOT.rglob("*")
        if path.is_file() and path.suffix.lower() in {".rdb", ".rdx"}
    ]
    assert forbidden == []


def test_managed_python_setup_files_are_present() -> None:
    assert (ROOT / "inst" / "python" / "provision_environment.R").is_file()
    assert (ROOT / "inst" / "python" / "managed-requirements.txt").is_file()
    assert (ROOT / "man" / "setup_importABCatlas.Rd").is_file()
    assert (ROOT / "TEST_AUTOMATIC_SETUP_9007.R").is_file()
    assert (ROOT / "TEST_FAILED_DATASETS_9007.R").is_file()


def test_readme_uses_automatic_python_setup() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "setup_importABCatlas()" in readme
    assert "You do **not** need to find a Python path" in readme
    assert "Sys.setenv(IMPORTABCATLAS_PYTHON" not in readme


def test_description_declares_reticulate_managed_runtime() -> None:
    description = (ROOT / "DESCRIPTION").read_text(encoding="utf-8")
    assert "reticulate (>= 1.41.0)" in description
    requirements = (ROOT / "inst" / "python" / "managed-requirements.txt").read_text(encoding="utf-8")
    assert "abc_atlas_access @ https://github.com/AllenInstitute/abc_atlas_access/archive/" in requirements
