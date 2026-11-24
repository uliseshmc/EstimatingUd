"""The coordinate-defined classes (utr5, utr3, introns) are ordinary non-cds regions."""

from pathlib import Path

import pytest

from estimating_ud import apps, regions

OLD_NAMES = ("introns5UTR", "introns3UTR", "introns_nonUTR")


def test_genic_regions_are_registered_as_noncds_regions():
    assert regions.GENIC_REGIONS == ["utr5", "utr3", "introns"]
    for name in regions.GENIC_REGIONS:
        assert name in regions.REGIONS
        assert name in regions.NONCDS_REGIONS
        assert name in regions.CONSTRAINT_REGIONS
        assert regions.raw_dir(Path("/data"), name, "22") == Path("/data") / name / "alldata_chrm22"
    assert len(regions.REGIONS) == 9 and len(set(regions.REGIONS)) == 9


def test_old_cut_at_mask_classes_are_gone():
    assert not set(OLD_NAMES) & set(regions.REGIONS)
    assert not hasattr(regions, "INTRON_REGIONS")
    for name in ("sample_UTR5", "sample_UTR3", "removeUTRs_fromintrons"):
        assert not hasattr(apps, name)


def test_coordinate_named_fragments_are_accepted(tmp_path):
    folder = tmp_path / "utr5" / "alldata_chrm22"
    folder.mkdir(parents=True)
    for name in ("homo_sapiens-22-100-150.fa", "homo_sapiens-22-18605814-18606514-1.fa", "homo_sapiens-X-5-9.fa"):
        (folder / name).write_text(">x\nA\n")
    regions.check_coordinate_named(folder, "utr5")


def test_gene_named_fragments_are_refused_for_genic_regions_only(tmp_path):
    folder = tmp_path / "raw"
    folder.mkdir()
    (folder / "homo_sapiens-22-100-150.fa").write_text(">x\nA\n")
    (folder / "ENSG00000100345.fa").write_text(">x\nA\n")
    for name in regions.GENIC_REGIONS:
        with pytest.raises(ValueError, match="ENSG00000100345.fa"):
            regions.check_coordinate_named(folder, name)
    regions.check_coordinate_named(folder, "intronsAR")  # per-gene names are right for the other regions
    regions.check_coordinate_named(tmp_path / "missing", "introns")  # nothing to refuse
