"""Check oligo orientation and indel handling independently of the report UI."""
import shutil
from pathlib import Path

import pytest

import primer_analysis as engine

PRIMER = "ACGTTGCAAGCTTAGCGATCGATGCTAGCAGTTCGACCTA"
FIXTURES = Path(__file__).resolve().parents[1] / "fixtures/reverse_primer"


def test_reverse_complement_preserves_iupac_meaning():
    assert engine.reverse_complement("ACGTURYSWKMBDHVN") == "NBDHVKMWSRYAACGT"
    assert engine.reverse_complement(PRIMER) == "TAGGTCGAACTGCTAGCATCGATCGCTAAGCTTGCAACGT"


@pytest.mark.parametrize("reverse", [False, True])
def test_missing_sequence_at_both_ends_is_padded_in_primer_orientation(reverse):
    middle = PRIMER[2:-2]
    # Fixture conversion is independent of the production reverse-complement helper.
    subject = middle.translate(str.maketrans("ACGT", "TGCA"))[::-1] if reverse else middle
    query, observed = engine.reconstruct_full_primer_alignment(
        PRIMER, middle, middle, 3, 38,
        len(middle) if reverse else 1, 1 if reverse else len(middle), subject,
    )
    assert query == PRIMER
    assert observed == "--" + middle + "--"
    assert engine.count_mismatches(query, observed) == 4
    assert engine.get_mismatch_positions(query, observed, 1, 40, PRIMER) == "1,2,39,40"


def test_reconstruction_keeps_query_gaps_even_without_subject_flanks():
    query, observed = engine.reconstruct_full_primer_alignment(
        "ACGTAC", "CG-TA", "CGATA", 2, 5, 2, 6, "",
    )
    assert query == "ACG-TAC"
    assert observed == "-CGATA-"
    assert engine.get_mismatch_details(query, observed, 1, 6, "ACGTAC") == "1:A>-,3:->A,6:C>-"


def test_indels_do_not_shift_positions_or_duplicate_chart_anchors():
    query, subject = "AC--GT", "ATACGA"
    assert engine.count_mismatches(query, subject) == 4
    assert engine.get_mismatch_positions(query, subject, 1, 4, "ACGT") == "2,4"
    assert engine.get_mismatch_details(query, subject, 1, 4, "ACGT") == "2:C>TAC,4:T>A"


def test_alignment_length_mismatch_cannot_silently_drop_trailing_bases():
    with pytest.raises(ValueError, match="same number of columns"):
        engine.count_mismatches("ACGT", "ACGTA")


@pytest.mark.skipif(not shutil.which("blastn"), reason="BLAST+ integration requires blastn")
def test_real_blast_is_invariant_to_subject_orientation():
    hits = engine.run_blastn(
        "audit_RIGHT", PRIMER, str(FIXTURES / "sequences.fasta"),
        execution=engine.BlastExecution(strict_errors=True, quiet=True),
    )
    expected = {
        "perfect": (0, ""), "first": (1, "1:A>T"), "last": (1, "40:A>C"),
        "both": (2, "1:A>T,40:A>C"), "internal": (1, "16:C>A"),
        "insertion": (1, "18:->G"), "deletion": (1, "19:T>-"),
    }
    assert len(hits) == 14
    for name, (count, details) in expected.items():
        forward, reverse = hits[name + "_plus"], hits[name + "_minus"]
        for hit in (forward, reverse):
            assert hit["mismatches"] == count, name
            assert hit["mismatch_details"] == details, name
            assert hit["pident"] == 100 * (40 - count) / 40, name
            assert hit["query_alignment"].replace("-", "") == PRIMER
            assert len(hit["query_alignment"]) == len(hit["subject_alignment"])
        assert forward["sstart"] < forward["send"]
        assert reverse["sstart"] > reverse["send"]
        for key in ("query_alignment", "subject_alignment", "mismatch_positions", "mismatch_details"):
            assert forward[key] == reverse[key], (name, key)
    # The full-length display includes terminal mismatches outside this local hit.
    assert (hits["both_minus"]["sstart"], hits["both_minus"]["send"]) == (79, 42)
    assert (hits["both_minus"]["qstart"], hits["both_minus"]["qend"]) == (2, 39)


@pytest.mark.skipif(not shutil.which("blastn"), reason="BLAST+ integration requires blastn")
def test_reverse_strand_iupac_and_terminal_mismatch_after_insertion(tmp_path):
    primer = PRIMER[:15] + "Y" + PRIMER[16:]
    # Y at position 16 matches C after reverse-complementing the subject.
    observed = PRIMER[:18] + "G" + PRIMER[18:-1] + "C"
    reverse = observed.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    fasta = tmp_path / "iupac.fasta"
    fasta.write_text(">sample\n" + "C" * 40 + reverse + "G" * 40 + "\n")
    hit = engine.run_blastn(
        "reverse", primer, str(fasta), execution=engine.BlastExecution(strict_errors=True, quiet=True),
    )["sample"]
    assert hit["sstart"] > hit["send"]
    assert hit["mismatch_details"] == "18:->G,40:A>C"
    assert hit["mismatches"] == 2
    assert hit["query_alignment"].replace("-", "") == primer
    assert hit["subject_alignment"] == observed
