from __future__ import annotations


def test_legacy_package_name_still_imports() -> None:
    import RNA_normalizer

    assert hasattr(RNA_normalizer, "PDBNormalizer")
    assert hasattr(RNA_normalizer, "PDBStruct")
