import json
import subprocess
import sys
from pathlib import Path

import pytest

import primer_checker
import primer_analysis
import primer_report


def test_iupac_matching_and_mismatch_counting():
    assert primer_checker.bases_match("Y", "T")
    assert primer_checker.bases_match("Y", "C")
    assert primer_checker.bases_match("R", "A")
    assert primer_checker.bases_match("N", "G")
    assert primer_checker.bases_match("U", "T")
    assert not primer_checker.bases_match("R", "C")
    assert primer_checker.count_mismatches("AY-", "ATG") == 1
    assert primer_checker.count_mismatches("YRNSWKMBDHV", "TGAGAGACACA") == 0


def test_mismatch_positions_include_partial_alignment_gaps():
    positions = primer_checker.get_mismatch_positions(
        qseq="CT",
        sseq="CA",
        qstart=2,
        qend=3,
        full_primer="ACTG",
    )
    assert positions == "1,3,4"
    details = primer_checker.get_mismatch_details(
        qseq="CT",
        sseq="CA",
        qstart=2,
        qend=3,
        full_primer="ACTG",
    )
    assert details == "1:A>-,3:T>A,4:G>-"


def test_iupac_matches_do_not_create_mismatch_positions_or_details():
    assert primer_checker.get_mismatch_positions("YRN", "TGA", 1, 3, "YRN") == ""
    assert primer_checker.get_mismatch_details("YRN", "TGA", 1, 3, "YRN") == ""


def test_reconstructed_alignment_shows_primer_end_bases_from_subject_sequence():
    primer = "AGCTACGATCGTACGATGCTACGTA"
    subject_bases = "AACTACGATCGTACGATGCTACGTA"
    subject = "N" * 235 + subject_bases + "N" * 10

    query_alignment, subject_alignment = primer_checker.reconstruct_full_primer_alignment(
        primer_seq=primer,
        qseq=primer[2:],
        sseq=subject_bases[2:],
        qstart=3,
        qend=len(primer),
        sstart=238,
        send=235 + len(primer),
        subject_sequence=subject,
    )

    assert query_alignment == primer
    assert subject_alignment == subject_bases
    assert primer_checker.get_mismatch_positions(query_alignment, subject_alignment, 1, len(primer), primer) == "2"
    assert primer_checker.get_mismatch_details(query_alignment, subject_alignment, 1, len(primer), primer) == "2:G>A"


def test_segment_parsing_supports_standard_and_contig_headers():
    assert primer_checker.get_segment("03-M|252500127") == "M"
    assert primer_checker.get_segment("01-HA|A/Netherlands/10685/2024") == "HA"
    assert primer_checker.get_segment("contig1|03-M|INFL16-2025") == "M"
    assert primer_checker.get_segment("Genome|252400716") == ""


def test_legacy_primer_validation_and_record_conversion():
    library = {
        "Influenza-B": {
            "triplex_InfB_F_NS": "CGTACCAYTTTTCTAGTAAAGT",
            "probe-VIC2_HA": "AGTACCATGGCAGATATGTCGHGGGTC",
        },
        "RSV-A": {"RSVQA1": "TTGGTGGGGCGACAGAAGCCAGGCGG"},
    }
    validation = primer_checker.validate_legacy_primer_library(library)
    assert validation.ok
    assert any("legacy format cannot store role metadata" in warning for warning in validation.warnings)

    records = primer_checker.legacy_library_to_records(library)
    assert records["Influenza-B"][0].segment == "NS"
    assert records["Influenza-B"][1].segment == "HA"
    assert records["Influenza-B"][1].role == "probe"


def test_legacy_primer_validation_rejects_invalid_sequences():
    validation = primer_checker.validate_legacy_primer_library({"RSV-A": {"bad": "ACGTZ"}})
    assert not validation.ok
    assert "unsupported IUPAC" in validation.errors[0]


def test_normalized_primer_validation_and_record_loading():
    normalized = {
        "schema_version": "1.0",
        "database_version": "2026-05-06",
        "schemes": [
            {
                "scheme_id": "dummy-influenza-b",
                "display_name": "Dummy Influenza B",
                "organism": "Influenza-B",
                "version": "2026-05-06",
                "status": "current",
                "source": "Dummy",
                "references": [],
                "primers": [
                    {
                        "id": "triplex_InfB_F_NS",
                        "name": "triplex_InfB_F_NS",
                        "sequence": "CGTACCAYTTTTCTAGTAAAGT",
                        "role": "forward_primer",
                        "segment": "NS",
                        "gene": "NS",
                        "pool": "triplex",
                        "strand": "plus",
                        "subtype_tags": [],
                        "notes": "",
                    }
                ],
            }
        ],
    }

    validation = primer_checker.validate_normalized_primer_library(normalized)
    assert validation.ok
    records = primer_checker.normalized_library_to_records(normalized)
    record = records["Influenza-B"][0]
    assert record.name == "triplex_InfB_F_NS"
    assert record.segment == "NS"
    assert record.role == "forward_primer"
    assert record.scheme_id == "dummy-influenza-b"
    assert record.database_version == "2026-05-06"


def test_load_primer_records_detects_normalized_json(tmp_path):
    primer_db = tmp_path / "normalized.json"
    primer_db.write_text(
        json.dumps(
            {
                "schema_version": "1.0",
                "database_version": "2026-05-06",
                "schemes": [
                    {
                        "scheme_id": "dummy-rsv-a",
                        "display_name": "Dummy RSV-A",
                        "organism": "RSV-A",
                        "version": "2026-05-06",
                        "status": "current",
                        "source": "Dummy",
                        "references": [],
                        "primers": [
                            {
                                "id": "RSVQA1",
                                "name": "RSVQA1",
                                "sequence": "TTGGTGGGGCGACAGAAGCCAGGCGG",
                                "role": "primer",
                                "segment": "",
                                "gene": "",
                                "pool": "",
                                "strand": "",
                                "subtype_tags": [],
                                "notes": "",
                            }
                        ],
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    records, validation = primer_checker.load_primer_records(str(primer_db))

    assert validation.ok
    assert records["RSV-A"][0].name == "RSVQA1"
    assert records["RSV-A"][0].scheme_id == "dummy-rsv-a"


def test_load_primer_records_supports_panel_bed_sequences(tmp_path):
    panel_dir = tmp_path / "panels" / "sars2-ngs-v1"
    panel_dir.mkdir(parents=True)
    (panel_dir / "amplicons.bed").write_text(
        "MN908947.3\t47\t78\tamp1_LEFT\t1\t+\tACGT\n"
        "MN908947.3\t419\t447\tamp1_RIGHT\t1\t-\tTGCA\n",
        encoding="utf-8",
    )
    primer_db = tmp_path / "panel_db.json"
    primer_db.write_text(
        json.dumps(
            {
                "schema_version": "2.0",
                "database_version": "2026-05-21",
                "panels": [
                    {
                        "panel_id": "sars2-ngs-v1",
                        "display_name": "SARS-CoV-2 NGS Panel v1",
                        "organism": "SARS-CoV-2",
                        "technology": "amplicon_ngs",
                        "panel_version": "1.0",
                        "source": "Dummy",
                        "reference": {"name": "MN908947.3", "coordinate_system": "0-based BED"},
                        "files": {"bed": "panels/sars2-ngs-v1/amplicons.bed"},
                        "mapping": {"pool_from_bed_field": 4},
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    records, validation = primer_checker.load_primer_records(str(primer_db))

    assert validation.ok
    by_name = {record.name: record for record in records["SARS-CoV-2"]}
    assert by_name["amp1_LEFT"].sequence == "ACGT"
    assert by_name["amp1_LEFT"].role == "forward_primer"
    assert by_name["amp1_RIGHT"].role == "reverse_primer"
    assert by_name["amp1_RIGHT"].pool == "1"
    assert by_name["amp1_RIGHT"].scheme_id == "sars2-ngs-v1"
    assert by_name["amp1_RIGHT"].reference_name == "MN908947.3"


def test_load_primer_records_supports_panel_fasta_with_amplicon_bed_mapping(tmp_path):
    panel_dir = tmp_path / "panels" / "sars2-ngs-v1"
    panel_dir.mkdir(parents=True)
    (panel_dir / "amplicons.bed").write_text("MN908947.3\t47\t447\tamp1\t2\t+\n", encoding="utf-8")
    (panel_dir / "primers.fasta").write_text(
        ">amp1_LEFT\nACGT\n>amp1_RIGHT\nTGCA\n",
        encoding="utf-8",
    )
    primer_db = tmp_path / "panel_db.json"
    primer_db.write_text(
        json.dumps(
            {
                "schema_version": "2.0",
                "database_version": "2026-05-21",
                "panels": [
                    {
                        "panel_id": "sars2-ngs-v1",
                        "display_name": "SARS-CoV-2 NGS Panel v1",
                        "organism": "SARS-CoV-2",
                        "technology": "amplicon_ngs",
                        "panel_version": "1.0",
                        "source": "Dummy",
                        "reference": {"name": "MN908947.3", "coordinate_system": "0-based BED"},
                        "files": {
                            "bed": "panels/sars2-ngs-v1/amplicons.bed",
                            "primers_fasta": "panels/sars2-ngs-v1/primers.fasta",
                        },
                        "mapping": {
                            "bed_name_field": "amplicon_id",
                            "primer_name_pattern": "{amplicon_id}_{role}",
                            "pool_from_bed_field": 4,
                        },
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    records, validation = primer_checker.load_primer_records(str(primer_db))

    assert validation.ok
    by_name = {record.name: record for record in records["SARS-CoV-2"]}
    assert set(by_name) == {"amp1_LEFT", "amp1_RIGHT"}
    assert by_name["amp1_LEFT"].pool == "2"
    assert by_name["amp1_LEFT"].role == "forward_primer"
    assert by_name["amp1_RIGHT"].role == "reverse_primer"
    assert by_name["amp1_RIGHT"].source_file == "panels/sars2-ngs-v1/primers.fasta"


def test_load_primer_records_supports_mixed_pcr_schemes_and_ngs_panels(tmp_path):
    panel_dir = tmp_path / "panels" / "sars2-ngs-v1"
    panel_dir.mkdir(parents=True)
    (panel_dir / "primers.bed").write_text(
        "MN908947.3\t47\t78\tngs_LEFT\t1\t+\tACGT\n",
        encoding="utf-8",
    )
    primer_db = tmp_path / "mixed_db.json"
    primer_db.write_text(
        json.dumps(
            {
                "schema_version": "2.0",
                "database_version": "2026-05-21",
                "schemes": [
                    {
                        "scheme_id": "dummy-sars-cov-2-pcr",
                        "display_name": "Dummy SARS-CoV-2 PCR",
                        "organism": "SARS-CoV-2",
                        "version": "2026-05-21",
                        "primers": [
                            {
                                "id": "triplex_SC2_F",
                                "name": "triplex_SC2_F",
                                "sequence": "CATATGAATTGAGGCAGAACGGT",
                                "role": "forward_primer",
                                "segment": "",
                            }
                        ],
                    }
                ],
                "panels": [
                    {
                        "panel_id": "sars2-ngs-v1",
                        "display_name": "SARS-CoV-2 NGS Panel v1",
                        "organism": "SARS-CoV-2",
                        "panel_version": "1.0",
                        "files": {"bed": "panels/sars2-ngs-v1/primers.bed"},
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    records, validation = primer_checker.load_primer_records(str(primer_db))

    assert validation.ok
    by_name = {record.name: record for record in records["SARS-CoV-2"]}
    assert by_name["triplex_SC2_F"].scheme_id == "dummy-sars-cov-2-pcr"
    assert by_name["ngs_LEFT"].scheme_id == "sars2-ngs-v1"
    assert by_name["ngs_LEFT"].pool == "1"


def test_load_primer_records_supports_virus_organized_pcr_and_ngs(tmp_path):
    panel_dir = tmp_path / "assets" / "SARS-CoV-2" / "Panel1"
    panel_dir.mkdir(parents=True)
    (panel_dir / "primers.bed").write_text(
        "MN908947.3\t47\t78\tngs_LEFT\t1\t+\tACGT\n",
        encoding="utf-8",
    )
    primer_db = tmp_path / "organized.json"
    primer_db.write_text(
        json.dumps(
            {
                "schema_version": "3.0",
                "database_version": "2026-05-21",
                "viruses": [
                    {
                        "organism": "SARS-CoV-2",
                        "pcr": {
                            "schemes": [
                                {
                                    "scheme_id": "dummy-sars-cov-2-pcr",
                                    "display_name": "Dummy SARS-CoV-2 PCR",
                                    "organism": "SARS-CoV-2",
                                    "version": "2026-05-21",
                                    "primers": [
                                        {
                                            "id": "triplex_SC2_F",
                                            "name": "triplex_SC2_F",
                                            "sequence": "CATATGAATTGAGGCAGAACGGT",
                                            "role": "forward_primer",
                                            "segment": "",
                                        }
                                    ],
                                }
                            ]
                        },
                        "ngs": {
                            "panels": [
                                {
                                    "panel_id": "sars2-ngs-v1",
                                    "display_name": "SARS-CoV-2 NGS Panel v1",
                                    "organism": "SARS-CoV-2",
                                    "panel_version": "1.0",
                                    "files": {"bed": "assets/SARS-CoV-2/Panel1/primers.bed"},
                                }
                            ]
                        },
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    records, validation = primer_checker.load_primer_records(str(primer_db))

    assert validation.ok
    by_name = {record.name: record for record in records["SARS-CoV-2"]}
    assert by_name["triplex_SC2_F"].scheme_id == "dummy-sars-cov-2-pcr"
    assert by_name["ngs_LEFT"].scheme_id == "sars2-ngs-v1"
    assert by_name["ngs_LEFT"].pool == "1"


def test_dummy_database_loads_unique_pcr_ngs_and_influenza_records():
    repo_root = Path(__file__).resolve().parents[1]
    records, validation = primer_checker.load_primer_records(str(repo_root / "primer_db/dummy_primers.json"))
    assert validation.ok
    assert not validation.warnings
    assert {organism: len(rows) for organism, rows in records.items()} == {"Demo-virus": 5, "Influenza-A": 8}
    all_records = [p for group in records.values() for p in group]
    assert len({(p.organism, p.scheme_id, p.name) for p in all_records}) == len(all_records)
    assert all(p.name.startswith("DUMMY_") and p.assay_name.startswith("DUMMY") for p in all_records)
    assert {p.database_version for p in all_records} == {"dummy-1.0"}
    for kind, count in [("pcr", 3), ("ngs", 2)]:
        virus, selected = primer_checker.select_primer_records(records, "Demo-virus", None, assay_type=kind)
        assert virus == "Demo-virus"
        assert len(selected) == count
        assert {p.assay_type for p in selected} == {kind}
    assert {p.segment for p in records["Influenza-A"]} == {"PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"}


def test_dummy_panel_database_loads_bed_with_separate_primer_fasta():
    repo_root = Path(__file__).resolve().parents[1]
    records, validation = primer_checker.load_primer_records(str(repo_root / "primer_db/assets/panel_database.example.json"))
    assert validation.ok
    assert len(records["Demo-virus"]) == 2
    by_name = {p.name: p for p in records["Demo-virus"]}
    left = by_name["DUMMY_LEFT"]
    assert left.sequence == "GTCAGACATCGATGCTACGTCAGGATCGTACCTAGCTGAC"
    assert left.pool == "1"
    assert left.start == "0"
    assert left.end == str(len(left.sequence))
    assert left.source_file == "dummy/primers.fasta"
    assert by_name["DUMMY_RIGHT"].role == "reverse_primer"


def test_legacy_to_normalized_conversion_preserves_names_and_sequences():
    legacy = {
        "Influenza-B": {
            "triplex_InfB_F_NS": "CGTACCAYTTTTCTAGTAAAGT",
            "probe-VIC2_HA": "AGTACCATGGCAGATATGTCGHGGGTC",
        }
    }

    normalized = primer_checker.legacy_library_to_normalized_database(
        legacy,
        database_version="2026-05-06",
        source="Dummy",
    )

    assert normalized["schema_version"] == "1.0"
    assert normalized["database_version"] == "2026-05-06"
    validation = primer_checker.validate_normalized_primer_library(normalized)
    assert validation.ok
    primers = normalized["schemes"][0]["primers"]
    assert {primer["name"] for primer in primers} == set(legacy["Influenza-B"])
    assert primers[0]["segment"] == "NS"
    assert primers[1]["role"] == "probe"


def test_influenza_subtype_selection_includes_h1_and_m_but_not_h3():
    records = primer_checker.legacy_library_to_records(
        {
            "Influenza-A": {
                "triplex_INFA_F1_M": "GGCTTGGAAACYTCCTTCTACCCGC",
                "H3_F_HA": "TCCTTTCAGYACGAGATAATC",
                "FluSw-H1-F236_HA": "TTCCGGTTATACGCTAGACTGAAGC",
            }
        }
    )
    selected = primer_checker.build_influenza_subtype_records(records, "H1")
    names = {record.name for record in selected}
    assert names == {"triplex_INFA_F1_M", "FluSw-H1-F236_HA"}


def test_influenza_subtype_selection_includes_matching_and_untagged_records():
    records = {
        "Influenza-A": [
            primer_checker.PrimerRecord(
                organism="Influenza-A",
                name="H3_F_HA",
                sequence="ACGT",
                segment="HA",
                subtype_tags=("H3",),
            ),
            primer_checker.PrimerRecord(
                organism="Influenza-A",
                name="H1_F_HA",
                sequence="ACGT",
                segment="HA",
                subtype_tags=("H1",),
            ),
            primer_checker.PrimerRecord(
                organism="Influenza-A",
                name="universal_M",
                sequence="ACGT",
                segment="M",
            ),
            primer_checker.PrimerRecord(
                organism="Influenza-A",
                name="untagged_HA",
                sequence="ACGT",
                segment="HA",
            ),
        ]
    }

    h3_names = {record.name for record in primer_checker.build_influenza_subtype_records(records, "H3")}
    h1_names = {record.name for record in primer_checker.build_influenza_subtype_records(records, "H1")}

    assert h3_names == {"H3_F_HA", "universal_M", "untagged_HA"}
    assert h1_names == {"H1_F_HA", "universal_M", "untagged_HA"}


def test_legacy_influenza_subtype_helper_uses_matching_and_untagged_names():
    legacy = {
        "Influenza-A": {
            "H3_F_HA": "ACGT",
            "H1_F_HA": "ACGT",
            "universal_M": "ACGT",
            "untagged_HA": "ACGT",
        }
    }

    assert set(primer_checker.build_influenza_subtype_primers(legacy, "H3")) == {
        "H3_F_HA",
        "universal_M",
        "untagged_HA",
    }


def test_filename_target_inference_for_batch_wrapper_names():
    assert primer_checker.infer_analysis_target_from_filename("H3.fasta") == primer_checker.AnalysisTarget(
        "influenza", "H3", "filename contains H3/H3N2"
    )
    assert primer_checker.infer_analysis_target_from_filename("run_H1N1_samples.fa").flu_type == "H1"
    assert primer_checker.infer_analysis_target_from_filename("INFB_test.fasta").flu_type == "B"
    assert primer_checker.infer_analysis_target_from_filename("SC2.fasta").virus_type == "SARS-CoV-2"
    assert primer_checker.infer_analysis_target_from_filename("RSVA_primercheck.fasta").virus_type == "RSV-A"
    assert primer_checker.infer_analysis_target_from_filename("RSVB.fasta").virus_type == "RSV-B"
    assert primer_checker.infer_analysis_target_from_filename(
        "Norovirus-GII_run1.fasta",
        available_organisms=["Norovirus-GII"],
    ).virus_type == "Norovirus-GII"
    assert primer_checker.infer_analysis_target_from_filename("unknown.fasta") is None


def test_batch_wrapper_help_shows_analysis_only_and_default_html_report():
    result = subprocess.run(
        [
            sys.executable,
            str(Path(__file__).resolve().parents[1] / "scripts" / "run_primer_checker_batch.py"),
            "--help",
        ],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert result.returncode == 0
    assert "--analysis-only" in result.stdout
    assert "batch_primer_report.html" in result.stdout


def test_batch_wrapper_dry_run_classifies_input_folder(tmp_path):
    input_folder = tmp_path / "fastas"
    input_folder.mkdir()
    for filename in ["H3.fasta", "SC2.fasta", "RSVB.fasta"]:
        (input_folder / filename).write_text(">sample\nACGT\n", encoding="utf-8")

    primer_db = tmp_path / "primers.json"
    primer_db.write_text(json.dumps({"SARS-CoV-2": {"primer-a": "ACGT"}, "Influenza-A": {"H3_HA": "ACGT"}}), encoding="utf-8")

    result = subprocess.run(
        [
            sys.executable,
            str(Path(__file__).resolve().parents[1] / "scripts" / "run_primer_checker_batch.py"),
            "--input-folder",
            str(input_folder),
            "--primers",
            str(primer_db),
            "--dry-run",
        ],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert result.returncode == 0
    assert "H3.fasta: influenza --flu-type H3" in result.stdout
    assert "SC2.fasta: SARS-CoV-2" in result.stdout
    assert "RSVB.fasta: RSV-B" in result.stdout


def test_process_fasta_filters_influenza_b_ns_records(monkeypatch, tmp_path):
    fasta = tmp_path / "influenza_b.fasta"
    fasta.write_text(
        ">01-HA|sample1\nACGT\n"
        ">08-NS|sample2\nACGT\n"
        ">03-M|sample3\nACGT\n",
        encoding="utf-8",
    )
    primer = primer_checker.PrimerRecord(
        organism="Influenza-B",
        name="triplex_InfB_F_NS",
        sequence="CGTACCAYTTTTCTAGTAAAGT",
        segment="NS",
    )
    monkeypatch.setattr(primer_analysis, "run_blastn", lambda *_args, **_kwargs: {})

    rows = primer_checker.process_fasta_file(str(fasta), "Influenza-B", [primer])

    assert len(rows) == 1
    assert rows[0]["Subject_Sequence_ID"] == "08-NS|sample2"
    assert rows[0]["Subject_Segment"] == "NS"
    assert rows[0]["Hit_Status"] == "no_hit"


def test_metadata_sample_matching_uses_exact_tokens_not_substrings():
    assert primer_checker.metadata_sample_matches_subject("454511", "Genome|454511")
    assert primer_checker.metadata_sample_matches_subject("Genome|4545", "Genome|4545")
    assert not primer_checker.metadata_sample_matches_subject("4545", "Genome|454511")
    assert not primer_checker.metadata_sample_matches_subject("Genome|4546", "Genome|4545")


def test_metadata_csv_enriches_result_rows(monkeypatch, tmp_path):
    fasta = tmp_path / "samples.fasta"
    fasta.write_text(">Genome|454511\nACGT\n>Genome|4545\nACGT\n", encoding="utf-8")
    metadata_csv = tmp_path / "metadata.csv"
    metadata_csv.write_text(
        "SampleID,Sample_Date,Ct_Value\n"
        "454511,2025-01-01,21.5\n"
        "Genome|4545,2025-01-08,24.0\n",
        encoding="utf-8",
    )
    records, validation = primer_checker.load_metadata_csv(str(metadata_csv))
    assert validation.ok

    primer = primer_checker.PrimerRecord(
        organism="RSV-A",
        name="primer-a",
        sequence="ACGT",
    )
    monkeypatch.setattr(primer_analysis, "run_blastn", lambda *_args, **_kwargs: {})

    rows = primer_checker.process_fasta_file(str(fasta), "RSV-A", [primer], metadata_records=records)

    by_subject = {row["Subject_Sequence_ID"]: row for row in rows}
    assert by_subject["Genome|454511"]["Metadata_Sample_ID"] == "454511"
    assert by_subject["Genome|454511"]["Sample_Date"] == "2025-01-01"
    assert by_subject["Genome|454511"]["Ct_Value"] == "21.5"
    assert by_subject["Genome|454511"]["Ct_Source"] == "Ct_Value"
    assert by_subject["Genome|4545"]["Metadata_Sample_ID"] == "Genome|4545"


def test_metadata_assay_specific_ct_columns_are_selected_by_primer_context(monkeypatch, tmp_path):
    fasta = tmp_path / "H3.fasta"
    fasta.write_text(">01-HA|454511\nACGT\n", encoding="utf-8")
    metadata_csv = tmp_path / "metadata.csv"
    metadata_csv.write_text(
        "SampleID,Sample_Date,CT_H1,CT_H3,Triplex-SC2_CT,Triplex-InfA_CT,Triplex-InfB_CT\n"
        "454511,2025-01-01,30.1,22.2,18.4,24.6,25.7\n",
        encoding="utf-8",
    )
    records, validation = primer_checker.load_metadata_csv(str(metadata_csv))
    assert validation.ok

    h3_primer = primer_checker.PrimerRecord(
        organism="Influenza-A",
        name="H3_F_HA",
        sequence="ACGT",
        segment="HA",
        subtype_tags=("H3",),
    )
    sc2_primer = primer_checker.PrimerRecord(
        organism="SARS-CoV-2",
        name="triplex_SC2_F",
        sequence="ACGT",
    )
    triplex_infa_primer = primer_checker.PrimerRecord(
        organism="Influenza-A",
        name="triplex_INFA_F1_M",
        sequence="ACGT",
        segment="M",
    )
    triplex_infb_primer = primer_checker.PrimerRecord(
        organism="Influenza-B",
        name="triplex_InfB_F_NS",
        sequence="ACGT",
        segment="NS",
    )
    monkeypatch.setattr(primer_analysis, "run_blastn", lambda *_args, **_kwargs: {})

    h3_rows = primer_checker.process_fasta_file(str(fasta), "Influenza-H3", [h3_primer], metadata_records=records)
    sc2_rows = primer_checker.process_fasta_file(str(fasta), "SARS-CoV-2", [sc2_primer], metadata_records=records)

    assert h3_rows[0]["Ct_Value"] == "22.2"
    assert h3_rows[0]["Ct_Source"] == "CT_H3"
    assert sc2_rows[0]["Ct_Value"] == "18.4"
    assert sc2_rows[0]["Ct_Source"] == "Triplex-SC2_CT"
    assert primer_checker.select_ct_value(records[0], "Influenza-H3", triplex_infa_primer) == ("24.6", "Triplex-InfA_CT")
    assert primer_checker.select_ct_value(records[0], "Influenza-B", triplex_infb_primer) == ("25.7", "Triplex-InfB_CT")


def test_html_report_contains_embedded_filterable_data_and_escapes_values():
    html = primer_report.build_html_report(
        [
            {
                "Fasta_File": "example.fasta",
                "Virus_Type": "Influenza-B",
                "Primer_Name": "triplex_InfB_F_NS",
                "Primer_Sequence": "CGTACCAYTTTTCTAGTAAAGT",
                "Primer_Segment": "NS",
                "Subject_Sequence_ID": "sample</script><b>",
                "Subject_Segment": "NS",
                "Hit_Status": "hit",
                "Percent_Identity": 95.5,
                "Alignment_Length": 22,
                "Mismatches": 1,
                "Gap_Openings": 0,
                "Query_Start": 1,
                "Query_End": 22,
                "Subject_Start": 10,
                "Subject_End": 31,
                "E_value": 0.001,
                "Bitscore": 40,
                "Mismatch_Positions": "7",
                "Mismatch_Details": "7:A>G",
                "Query_Alignment": "CGTACCAYTTTTCTAGTAAAGT",
                "Subject_Alignment": "TCCTCAGYTCACTCTTCGAGCG",
                "Metadata_Sample_ID": "sample-1",
                "Sample_Date": "2025-01-01",
                "Ct_Value": "22.4",
                "Ct_Source": "CT_H3",
            }
        ],
        previous_reports=[
            {
                "name": "old_report.csv",
                "path": "/tmp/old_report.csv",
                "warnings": [],
                "rows": [
                    {
                        "Fasta_File": "old.fasta",
                        "Virus_Type": "Influenza-B",
                        "Primer_Name": "triplex_InfB_F_NS",
                        "Primer_Sequence": "CGTACCAYTTTTCTAGTAAAGT",
                        "Primer_Segment": "NS",
                        "Subject_Sequence_ID": "old-sample",
                        "Subject_Segment": "NS",
                        "Hit_Status": "hit",
                        "Percent_Identity": 90,
                        "Mismatches": 2,
                        "Mismatch_Positions": "1,2",
                        "Mismatch_Details": "1:T>C,2:C>T",
                    }
                ],
            }
        ],
    )
    english_text = primer_report.load_report_translations()["en"]

    assert "filter-primer" in html
    assert "filter-risk" in html
    assert english_text["risk"] in html
    assert "filter-status" not in html
    assert "filter-mismatches" not in html
    assert "Maximum mismatches" not in html
    assert english_text["primer_investigation_title"] in html
    assert english_text["ngs_title"] in html
    assert "renderNgsPanelOverview" in html
    assert english_text["no_viable_primer"] in html
    assert english_text["ngs_note"] in html
    assert "Review level" in html
    assert "RISK_THRESHOLDS" in html
    assert english_text["mismatch_distribution_title"] in html
    assert "<h2>Mismatch Count Distribution</h2>" not in html
    assert english_text["chart_title"] in html
    assert "Mismatch distribution along primer sequence" not in html
    assert "Mismatch_Details" in html
    assert "A&gt;G" in html or "A\\u003eG" in html
    assert "percent_axis" in html
    assert "sampleBaseFromChange" in html
    assert english_text["sample_base"] in html
    assert "data-language=\"no\"" in html
    assert "Norsk" in html
    assert "Sample timeline" not in html
    assert english_text["modal_title"] in html
    assert "alignment-modal" in html
    assert "data-alignment-key" in html
    assert english_text["no_alignment"] in html
    assert "IUPAC_BASES_FOR_REPORT" in html
    assert "Y: ['C', 'T']" in html
    assert "if (!state.sortKey)" in html
    assert english_text["doc_title"] in html
    assert 'data-tab="documentation"' in html
    assert 'class="doc-lead"' in html
    assert 'class="doc-faq-list"' in html
    assert 'class="doc-risk-list"' in html
    assert 'class="doc-two-column"' in html
    assert 'class="doc-step-list"' in html
    assert 'class="report-note"' in html
    assert english_text["doc_faq_1_q"] in html
    assert english_text["doc_workflow_title"] in html
    assert english_text["overview_note_label"] in html
    assert english_text["chart_note_what_label"] in html
    assert english_text["col_Terminal_Mismatch_Share"] in html
    assert english_text["two_plus_rate"] in html
    assert english_text["three_plus_rate"] in html
    assert "t('two_plus_rate') + '<strong><br>' + formatPercent(stats.twoPlusMismatchRate)" in html
    assert "t('three_plus_rate') + '<strong><br>' + formatPercent(stats.threePlusMismatchRate)" in html
    assert english_text["doc_faq_5"] in html
    assert english_text["doc_risk_high"] in html
    assert "terminalMismatchHits / hitRows.length" in html
    assert "highTerminalMismatchRate: 0.05" in html
    assert english_text["sample_rows_sorted"] in html
    assert "doc_faq_6_q" in html
    assert "formatIdentityValue" in html
    assert "toFixed(2)" in html
    assert "Ct_Value" not in html
    assert "Ct_Source" not in html
    assert "Sample_Date" not in html
    assert "Metadata_Sample_ID" not in html
    assert "2025-01-01" not in html
    assert "Avg_Mismatches" not in html
    assert "Average mismatches" not in html
    assert "Primers to Review" not in html
    assert "renderReviewPrimerList" not in html
    assert "data-detail-sort" in html
    assert "data-table-toggle" in html
    assert "primerPanelId" in html
    assert "primer-link" in html
    assert english_text["jump_to_primer"] in html
    assert "Showing top 10" not in html
    assert english_text["show_top_10"] in html
    assert "previous-report-data" not in html
    assert "Attached Previous Reports" not in html
    assert "old_report.csv" not in html
    assert "triplex_InfB_F_NS" in html
    assert "sample</script>" not in html
    assert "\\u003c/script\\u003e" in html


def test_report_text_files_are_loaded_and_validated(tmp_path):
    english_path = tmp_path / "english.json"
    norwegian_path = tmp_path / "norwegian.json"
    english = {
        "report_title": "Custom report",
        "open_alignment": "Open alignment for {sample}",
    }
    norwegian = {
        "report_title": "Tilpasset rapport",
        "open_alignment": "Åpne sekvenssammenstilling for {sample}",
    }
    english_path.write_text(json.dumps(english), encoding="utf-8")
    norwegian_path.write_text(json.dumps(norwegian), encoding="utf-8")

    translations = primer_report.load_report_translations(tmp_path)

    assert translations["en"]["report_title"] == "Custom report"
    assert translations["no"]["report_title"] == "Tilpasset rapport"

    norwegian["open_alignment"] = "Åpne sekvenssammenstilling"
    norwegian_path.write_text(json.dumps(norwegian), encoding="utf-8")

    with pytest.raises(ValueError, match="Placeholders.*open_alignment"):
        primer_report.load_report_translations(tmp_path)


def test_previous_report_csv_loader_reads_rows_and_warns_on_bad_csv(tmp_path):
    previous_csv = tmp_path / "previous.csv"
    previous_csv.write_text(
        "Fasta_File,Virus_Type,Primer_Name,Mismatches,Hit_Status\n"
        "old.fasta,Influenza-B,primer-a,2,hit\n",
        encoding="utf-8",
    )

    report = primer_report.load_previous_report_csv(str(previous_csv))

    assert report["name"] == "previous.csv"
    assert report["rows"][0]["Primer_Name"] == "primer-a"
    assert report["rows"][0]["Mismatches"] == "2"

    missing_report = primer_report.load_previous_report_csv(str(tmp_path / "missing.csv"))
    assert missing_report["warnings"]


def test_html_report_build_failure_preserves_existing_output(monkeypatch, tmp_path):
    output = tmp_path / "report.html"
    output.write_text("existing report", encoding="utf-8")

    def fail_to_build(*_args, **_kwargs):
        raise ValueError("invalid report text")

    monkeypatch.setattr(primer_report, "build_html_report", fail_to_build)
    primer_report.write_html_report([{"Primer_Name": "primer-a"}], str(output))

    assert output.read_text(encoding="utf-8") == "existing report"

def test_run_blastn_cleans_query_file_when_blast_is_missing(monkeypatch, tmp_path):
    query_file = tmp_path / "query.fasta"

    class FakeTempFile:
        name = str(query_file)

        def __enter__(self):
            self.handle = query_file.open("w", encoding="utf-8")
            return self.handle

        def __exit__(self, exc_type, exc, tb):
            self.handle.close()

    monkeypatch.setattr(primer_checker.tempfile, "NamedTemporaryFile", lambda **_kwargs: FakeTempFile())

    def missing_blast(*_args, **_kwargs):
        raise FileNotFoundError("blastn")

    monkeypatch.setattr(primer_checker.subprocess, "run", missing_blast)

    with pytest.raises(SystemExit):
        primer_checker.run_blastn("primer", "ACGT", "subject.fasta")

    assert not query_file.exists()


def test_validate_primers_cli_does_not_require_virus_or_fasta(tmp_path):
    primer_db = tmp_path / "primers.json"
    primer_db.write_text(json.dumps({"RSV-A": {"RSVQA1": "ACGT"}}), encoding="utf-8")

    result = subprocess.run(
        [
            sys.executable,
            str(Path(__file__).resolve().parents[1] / "primer_checker.py"),
            "--primers",
            str(primer_db),
            "--validate-primers",
        ],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    assert result.returncode == 0
    assert "validation succeeded" in result.stdout
