"""Database-driven influenza selection must not widen to unrelated tagged primers."""

import pytest

import primer_analysis as engine


@pytest.mark.parametrize("tag", ["H2", "H5N1", "H7N9", "H9N2", "H10", "H16", "H18N11", "CUSTOM-LINEAGE"])
def test_arbitrary_subtype_selects_exact_tag_and_shared_primers(tag):
    records = {"Influenza-A": [
        engine.PrimerRecord("Influenza-A", "specific", "ACGT", subtype_tags=(tag,), segment="PB2"),
        engine.PrimerRecord("Influenza-A", "shared", "ACGT", segment="NA"),
        engine.PrimerRecord("Influenza-A", "other", "ACGT", subtype_tags=("H3",), segment="NP"),
    ]}
    assert set(engine.influenza_selections(records)) == {"A", tag, "H3"}
    virus, selected = engine.select_primer_records(records, "influenza", tag.lower())
    assert virus == f"Influenza-{tag}"
    assert {p.name for p in selected} == {"specific", "shared"}
    with pytest.raises(SystemExit):
        engine.select_primer_records(records, "influenza", "NOT-IN-DATABASE")


def test_types_and_identical_lineage_tags_never_mix_organisms():
    records = {
        organism: [engine.PrimerRecord(organism, organism, "ACGT", segment="HEF", subtype_tags=("LINEAGE1",))]
        for organism in ["Influenza-A", "Influenza-B", "Influenza-C", "Influenza-D"]
    }
    for selection, organism in [("LINEAGE1", "Influenza-A"), ("B/LINEAGE1", "Influenza-B"),
                                ("C/LINEAGE1", "Influenza-C"), ("D", "Influenza-D")]:
        _, selected = engine.select_primer_records(records, "influenza", selection)
        assert [p.organism for p in selected] == [organism]


def test_legacy_inference_supports_all_standard_segments_and_full_subtype_tokens():
    for segment in ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS", "HEF"]:
        records = engine.legacy_library_to_records({"Influenza-A": {f"H5N1_F_{segment}": "ACGT"}})
        primer = records["Influenza-A"][0]
        assert primer.segment == segment
        assert primer.subtype_tags == ("H5N1",)
        assert engine.filter_subject_ids_for_primer(
            [f"01-{segment}|right", "01-OTHER|wrong"], primer, "Influenza-H5N1", "test.fa"
        ) == [f"01-{segment}|right"]
    assert engine.infer_subtype_tags("H10_F_PB2") == ("H10",)


def test_exact_tags_do_not_expand_h5_into_h5n1_and_empty_tags_override_names():
    assert engine.read_subtype_tags({"subtype_tags": []}, "H3_NA") == ()
    assert engine.read_subtype_tags({}, "H5N1_PB2") == ("H5N1",)
    assert engine.read_subtype_tags({"subtype_tags": ["h5n1", "H5N1", "h7n9"]}, "H3_PB2") == ("H5N1", "H7N9")
    records = engine.legacy_library_to_records({"Influenza-A": {"H5_HA": "ACGT", "H5N1_NA": "ACGT"}})
    assert [p.name for p in engine.build_influenza_subtype_records(records, "H5")] == ["H5_HA"]


def test_batch_uses_database_tags_and_preserves_family_fallback():
    records = engine.legacy_library_to_records({"Influenza-A": {"H5N1_PB2": "ACGT", "H1_HA": "ACGT"}})
    assert engine.infer_analysis_target_from_filename("run_H5N1.fa", primer_records=records).flu_type == "H5N1"
    assert engine.infer_analysis_target_from_filename("run_H1N1.fa", primer_records=records).flu_type == "H1"
    assert engine.infer_analysis_target_from_filename("run_H3.fa", primer_records=records) is None
    assert engine.infer_analysis_target_from_filename("run_H1_H5N1.fa", primer_records=records) is None
