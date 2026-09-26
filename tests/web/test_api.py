import csv
import io
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

import primer_analysis as engine
from api.index import app
from web_service import service

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "fixtures/simple_test_data"


@pytest.fixture
def client(monkeypatch):
    monkeypatch.setenv("PRIMER_DATABASE_PATH", str(FIXTURES / "simple_primers.json"))
    with TestClient(app, raise_server_exceptions=False) as client:
        yield client


def post(client, data=None, content=None, name="SARS.fasta", metadata=None):
    files = [("files", (name, content if content is not None else (FIXTURES / "SARS.fasta").read_bytes()))]
    if metadata is not None:
        files.append(("metadata", ("metadata.csv", metadata)))
    return client.post("/api/analyze", data=data or {"virus": "SARS-CoV-2", "assay_type": "pcr"}, files=files)


def test_catalog_matches_real_database(client):
    response = client.get("/api/catalog")
    assert response.status_code == 200
    assert "influenza" in [v["id"] for v in response.json()["viruses"]]
    assert len(response.json()["database"]["sha256"]) == 64


@pytest.mark.parametrize(
    "content", [b"", b"ACGT", b">x\n", b">x\nACGTZ", b">\nACGT", b">x\nACGT\n>x\nACGT", b">x\nA-CG", b"\xff"]
)
def test_invalid_fasta(client, content):
    response = post(client, content=content)
    assert response.status_code == 422
    assert "Traceback" not in response.text


def test_missing_file(client):
    response = client.post("/api/analyze", data={"virus": "RSV-A"})
    assert response.status_code == 422


@pytest.mark.parametrize(
    "selection",
    [
        {"virus": "unknown"},
        {"virus": "influenza"},
        {"virus": "influenza", "flu_type": "H5"},
        {"virus": "RSV-A", "flu_type": "H3"},
        {"virus": "SARS-CoV-2", "assay_type": "invalid"},
        {"virus": "SARS-CoV-2", "assay_id": "unknown"},
    ],
)
def test_invalid_selection(client, selection):
    assert post(client, data=selection).status_code == 422


@pytest.mark.parametrize(
    "metadata", [b"A,B\n1,2", b"SampleID,Ct\nx,2,3", b'SampleID,Ct\nx,"oops', b"SampleID,SampleID\nx,x"]
)
def test_invalid_metadata(client, metadata):
    assert post(client, metadata=metadata).status_code == 422


def test_upload_and_response_limits(client, monkeypatch):
    monkeypatch.setattr(service, "MAX_REQUEST_BYTES", 100)
    assert post(client).status_code == 413


def test_chunked_body_limit(client, monkeypatch):
    monkeypatch.setattr(service, "MAX_REQUEST_BYTES", 100)
    response = client.post(
        "/api/analyze",
        content=iter([b"a" * 60, b"b" * 60]),
        headers={"content-type": "multipart/form-data; boundary=x"},
    )
    assert response.status_code == 413


def test_sanitized_filenames():
    assert service.safe_filename("../../sample.fasta", {".fasta"}) == "sample.fasta"
    assert service.safe_filename("C:\\private\\sample.fasta", {".fasta"}) == "sample.fasta"


def test_database_failure_hides_paths(client, monkeypatch):
    monkeypatch.setenv("PRIMER_DATABASE_PATH", "/private/secret/missing.json")
    response = client.get("/api/catalog")
    assert response.status_code == 503
    assert "/private" not in response.text


@pytest.mark.parametrize(
    "failure,status,code",
    [
        (engine.BlastError("secret"), 502, "blast_failure"),
        (engine.BlastUnavailableError("secret"), 503, "blast_unavailable"),
        (engine.AnalysisTimeoutError("Analysis timed out."), 504, "analysis_timeout"),
        (ValueError("/private/internal"), 500, "analysis_error"),
    ],
)
def test_safe_errors(client, monkeypatch, failure, status, code):
    def fail(*args, **kwargs):
        raise failure

    monkeypatch.setattr(service, "analyze", fail)
    response = post(client)
    assert response.status_code == status
    assert response.json()["error"]["code"] == code
    assert "secret" not in response.text and "/private" not in response.text


@pytest.mark.skipif(not shutil.which("blastn"), reason="BLAST+ integration requires blastn")
@pytest.mark.parametrize(
    "filename,virus,subtype",
    [("SARS.fasta", "SARS-CoV-2", None), ("H3.fasta", "influenza", "H3"), ("RSVA.fasta", "RSV-A", None)],
)
def test_real_blast_cli_web_equivalence(client, tmp_path, filename, virus, subtype, capsys):
    options = {"virus": virus, "assay_type": "all"}
    args = [
        sys.executable,
        "primer_checker.py",
        "--primers",
        str(FIXTURES / "simple_primers.json"),
        "--virus",
        virus,
        "--fasta",
        str(FIXTURES / filename),
        "--metadata-csv",
        str(FIXTURES / "metadata.csv"),
        "--output",
        str(tmp_path / "cli.csv"),
        "--html-report",
        str(tmp_path / "cli.html"),
    ]
    if subtype:
        options["flu_type"] = subtype
        args += ["--flu-type", subtype]
    subprocess.run(args, cwd=ROOT, check=True, capture_output=True)
    response = post(
        client, options, (FIXTURES / filename).read_bytes(), filename, (FIXTURES / "metadata.csv").read_bytes()
    )
    assert response.status_code == 200, response.text
    result = response.json()
    assert list(csv.DictReader(io.StringIO(result["downloads"]["csv"]))) == list(
        csv.DictReader((tmp_path / "cli.csv").open())
    )
    assert result["summary"]["hits"] == 1
    assert result["manifest"]["blast"]["word_size"] == 4
    assert "Analysis provenance" in result["downloads"]["html"]
    assert json.loads(result["downloads"]["manifest"]) == result["manifest"]
    assert response.headers["cache-control"] == "no-store"
    assert not capsys.readouterr().out


def test_blast_deadline_and_query_cleanup(tmp_path, monkeypatch):
    monkeypatch.setattr(engine.tempfile, "tempdir", str(tmp_path))

    def timed_out(*args, **kwargs):
        raise subprocess.TimeoutExpired("blastn", 1)

    monkeypatch.setattr(engine.subprocess, "run", timed_out)
    with pytest.raises(engine.AnalysisTimeoutError):
        engine.run_blastn(
            "p", "ACGT", "subject.fasta", execution=engine.BlastExecution("blastn", time.monotonic() + 1, True)
        )
    assert not list(tmp_path.iterdir())


def test_blast_failure_is_not_a_no_hit(monkeypatch):
    monkeypatch.setattr(
        engine.subprocess, "run", lambda *args, **kwargs: subprocess.CompletedProcess(args, 1, "", "private details")
    )
    with pytest.raises(engine.BlastError, match="BLAST analysis failed"):
        engine.run_blastn("p", "ACGT", "subject.fasta", execution=engine.BlastExecution("blastn", strict_errors=True))


def test_executable_resolution_precedence(tmp_path, monkeypatch):
    binary = tmp_path / "blastn"
    binary.write_text("#!/bin/sh\nexit 0\n")
    binary.chmod(0o755)
    monkeypatch.setenv("BLASTN_PATH", str(binary))
    assert engine.resolve_blastn() == str(binary)
    monkeypatch.setenv("BLASTN_PATH", str(tmp_path / "missing"))
    with pytest.raises(engine.BlastUnavailableError):
        engine.resolve_blastn()


@pytest.fixture
def simulated_blast(monkeypatch):
    """Exercise selection/validation without paying per-primer process startup."""
    monkeypatch.setattr(engine, "resolve_blastn", lambda: "/test/blastn")
    monkeypatch.setattr(
        engine.subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(a, 0, "blastn: test\n", "")
    )
    monkeypatch.setattr(engine, "run_blastn", lambda *a, **k: {})


@pytest.mark.parametrize("subtype", ["A", "H1", "H3", "B"])
def test_influenza_subtypes_use_canonical_selection(client, monkeypatch, simulated_blast, subtype):
    monkeypatch.delenv("PRIMER_DATABASE_PATH")
    response = post(
        client,
        {"virus": "influenza", "flu_type": subtype, "assay_type": "pcr"},
        b">01-HA|sample\nACGTACGT\n>03-M|sample\nACGTACGT\n>08-NS|sample\nACGTACGT",
    )
    assert response.status_code == 200, response.text
    db, _, _ = service.load_database()
    expected_virus, expected = engine.select_primer_records(db, "influenza", subtype, "pcr")
    assert response.json()["summary"]["primers"] == len(expected)
    assert {r["Virus_Type"] for r in response.json()["rows"]} == {expected_virus}


@pytest.mark.parametrize(
    "assay_type,assay_id,count",
    [
        ("pcr", "fhi-sars-cov-2", 3),
        ("ngs", "sars2-ngs-vmidt-2.2", 68),
        ("all", None, 268),
    ],
)
def test_pcr_ngs_and_all_selection(client, monkeypatch, simulated_blast, assay_type, assay_id, count):
    monkeypatch.delenv("PRIMER_DATABASE_PATH")
    selection = {"virus": "SARS-CoV-2", "assay_type": assay_type}
    if assay_id:
        selection["assay_id"] = assay_id
    response = post(client, selection)
    assert response.status_code == 200, response.text
    assert response.json()["summary"]["primers"] == count
    if assay_id:
        assert {r["Assay_ID"] for r in response.json()["rows"]} == {assay_id}


def test_multiple_files_and_cleanup_on_failure(client, monkeypatch, simulated_blast):
    observed = []

    def inspect_files(*args, **kwargs):
        observed.append(Path(args[2]))
        assert observed[-1].is_file()
        return {}

    monkeypatch.setattr(engine, "run_blastn", inspect_files)
    response = client.post(
        "/api/analyze",
        data={"virus": "SARS-CoV-2"},
        files=[("files", ("../../one.fasta", b">x\nACGT")), ("files", ("two.fasta", b">x\nACGT"))],
    )
    assert response.status_code == 200
    assert response.json()["summary"]["files"] == 2
    assert {r["Fasta_File"] for r in response.json()["rows"]} == {"one.fasta", "two.fasta"}
    assert all(not p.exists() and not p.parent.parent.exists() for p in observed)

    def fail(*args, **kwargs):
        observed.append(Path(args[2]))
        raise engine.BlastError("failure")

    monkeypatch.setattr(engine, "run_blastn", fail)
    assert post(client).status_code == 502
    assert all(not p.exists() and not p.parent.parent.exists() for p in observed)


def test_large_results_are_rejected_and_work_is_bounded(client, monkeypatch, simulated_blast):
    monkeypatch.setattr(service, "MAX_RESPONSE_BYTES", 50)
    assert post(client).json()["error"]["code"] == "result_too_large"
    monkeypatch.setattr(service, "MAX_COMPARISONS", 0)
    assert post(client).json()["error"]["code"] == "work_limit"


def test_downloads_escape_untrusted_content(client, simulated_blast):
    response = post(client, content=b">=1+1\nACGT\n></script><script>alert(1)</script>\nACGT")
    assert response.status_code == 200
    result = response.json()
    csv_rows = list(csv.DictReader(io.StringIO(result["downloads"]["csv"])))
    assert csv_rows[0]["Subject_Sequence_ID"] == "'=1+1"
    assert result["rows"][0]["Subject_Sequence_ID"] == "=1+1"
    assert "</script><script>alert(1)</script>" not in result["downloads"]["html"]


def test_missing_segments_produce_actionable_error(client, simulated_blast):
    response = post(client, {"virus": "influenza", "flu_type": "H3"})
    assert response.status_code == 422
    assert response.json()["error"]["code"] == "no_comparisons"


def test_trimmed_fasta_headers_match_both_engine_readers(client, simulated_blast):
    response = post(client, content=b"  >sample  \n ACGT \n")
    assert response.status_code == 200
    assert response.json()["rows"][0]["Subject_Sequence_ID"] == "sample"
