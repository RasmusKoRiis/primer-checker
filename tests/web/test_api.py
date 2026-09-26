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
    "template,fasta,selection,records,comparisons",
    [
        ("normalized", "general", {"virus": "Example-virus"}, 2, 4),
        ("legacy", "general", {"virus": "Example-virus"}, 2, 4),
        ("influenza", "influenza", {"virus": "influenza", "flu_type": "H1"}, 8, 2),
        ("influenza", "influenza", {"virus": "influenza", "flu_type": "H5N1"}, 8, 3),
        ("influenza", "influenza", {"virus": "influenza", "flu_type": "A"}, 8, 4),
    ],
)
def test_downloadable_format_templates(client, template, fasta, selection, records, comparisons):
    # The UI previews and downloads this same data, so examples exercise the real parser.
    examples = json.loads((ROOT / "app/lib/input-examples.json").read_text())
    database = json.dumps(examples["databases"][template]).encode()
    assert client.post("/api/database", files={"database": ("template.json", database)}).status_code == 200
    response = client.post(
        "/api/preflight",
        data={**selection, "assay_type": "pcr"},
        files=[
            ("database", ("template.json", database)),
            ("files", ("template.fasta", examples["fasta"][fasta].encode())),
        ],
    )
    assert response.status_code == 200, response.text
    assert response.json()["warnings"] == []
    assert response.json()["workload"]["records"] == records
    assert response.json()["workload"]["comparisons"] == comparisons


@pytest.mark.parametrize("header", ["sample_001 HA", "sample_001|HA"])
def test_documented_incorrect_influenza_headers_do_not_match(client, header):
    examples = json.loads((ROOT / "app/lib/input-examples.json").read_text())
    response = client.post(
        "/api/preflight",
        data={"virus": "influenza", "flu_type": "H1", "assay_type": "pcr"},
        files=[
            ("database", ("template.json", json.dumps(examples["databases"]["influenza"]).encode())),
            ("files", ("template.fasta", f">{header}\nACGTACGT\n".encode())),
        ],
    )
    assert response.status_code == 422
    assert response.json()["error"]["code"] == "no_comparisons"


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
    web_rows = list(csv.DictReader(io.StringIO(result["downloads"]["csv"])))
    assert [{k: r[k] for k in engine.CSV_FIELDNAMES} for r in web_rows] == list(
        csv.DictReader((tmp_path / "cli.csv").open())
    )
    assert web_rows[0]["Database_SHA256"] == result["manifest"]["database"]["sha256"]
    assert json.loads(web_rows[0]["BLAST_Config_JSON"])["word_size"] == 4
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


def custom_database(assay_type="pcr"):
    return {
        "schema_version": "1.0",
        "database_version": "test-1",
        "schemes": [
            {
                "scheme_id": "my-scheme",
                "display_name": "My primers",
                "organism": "Custom-virus",
                "version": "1",
                "assay_type": assay_type,
                "primers": [
                    {
                        "id": "target_F",
                        "name": "target_F",
                        "role": "forward",
                        "segment": "",
                        "sequence": "ACGTTGCAAGCTTAGCGATCGATGCTAGCA",
                        "pool": "1",
                    }
                ],
            }
        ],
    }


def custom_request(client, endpoint, database=None, fasta=b">sample\nACGTTGCAAGCTTAGCGATCGATGCTAGCA\n", **selection):
    return client.post(
        endpoint,
        data={"virus": "Custom-virus", **selection},
        files=[
            ("files", ("custom.fasta", fasta)),
            ("database", ("custom.json", json.dumps(database or custom_database()).encode())),
        ],
    )


@pytest.mark.parametrize("assay_type", ["pcr", "ngs"])
def test_uploaded_database_is_selected_and_fingerprinted(client, simulated_blast, tmp_path, assay_type):
    database = custom_database(assay_type)
    content = json.dumps(database).encode()
    uploaded = client.post("/api/database", files={"database": ("custom.json", content)})
    assert uploaded.status_code == 200, uploaded.text
    catalog = uploaded.json()
    assert catalog["assays"][0]["type"] == assay_type
    assert catalog["viruses"] == [{"id": "Custom-virus", "name": "Custom-virus", "subtypes": []}]
    response = custom_request(client, "/api/analyze", database, assay_type=assay_type)
    assert response.status_code == 200, response.text
    result = response.json()
    assert result["manifest"]["database"] == catalog["database"]
    assert result["manifest"]["database"]["source"] == "uploaded"
    assert result["summary"]["primers"] == 1
    assert result["rows"][0]["Assay_Type"] == assay_type
    # Downloaded builder schema must work in the original CLI loader as well.
    path = tmp_path / "custom.json"
    path.write_bytes(content)
    cli_records, _ = engine.load_primer_records(str(path))
    _, selected = engine.select_primer_records(cli_records, "Custom-virus", None, assay_type, "my-scheme")
    assert selected[0].sequence == result["rows"][0]["Primer_Sequence"]
    # Request-scoped database must not replace the installed catalog.
    assert "Custom-virus" not in [v["id"] for v in client.get("/api/catalog").json()["viruses"]]


@pytest.mark.parametrize(
    "content",
    [
        b"{}",
        b"[]",
        b"not-json",
        b'{"virus":{"p":"ACGT","p":"AAAA"}}',
        b'{"virus":{"p\\n>injected":"ACGT"}}',
        b'{"virus":{"p":123}}',
        b'{"schema_version":"1.0","panels":[{"files":{"primers_fasta":"/private/secrets"}}]}',
        b'{"schema_version":"1.0","viruses":[]}',
    ],
)
def test_malformed_custom_database_is_rejected_without_asset_reads(client, monkeypatch, content):
    monkeypatch.setattr(engine, "load_primer_records", lambda *_: pytest.fail("Uploaded database resolved a path"))
    response = client.post("/api/database", files={"database": ("custom.json", content)})
    assert response.status_code == 422, response.text
    assert "Traceback" not in response.text and "/private" not in response.text


def test_custom_database_limits_and_legacy(client):
    assert client.post("/api/database", files={"database": ("custom.json", b" " * 250_001)}).status_code == 413
    database = custom_database()
    database["schemes"][0]["primers"][0]["sequence"] = "A" * 201
    assert client.post("/api/database", files={"database": ("custom.json", json.dumps(database))}).status_code == 422
    database = custom_database()
    database["schemes"][0]["primers"] = [
        dict(database["schemes"][0]["primers"][0], id=str(i), name=str(i)) for i in range(501)
    ]
    assert client.post("/api/database", files={"database": ("custom.json", json.dumps(database))}).status_code == 422
    response = client.post(
        "/api/database", files={"database": ("legacy.json", json.dumps({"Custom-virus": {"target_F": "ACGT"}}))}
    )
    assert response.status_code == 200
    assert response.json()["assays"][0]["primers"] == 1


def test_preflight_counts_work_without_starting_blast(client, monkeypatch, simulated_blast):
    monkeypatch.setattr(service, "blast_version", lambda: pytest.fail("Preflight launched BLAST"))
    response = custom_request(client, "/api/preflight")
    assert response.status_code == 200, response.text
    assert response.json()["workload"] == {
        "files": 1,
        "records": 1,
        "primers": 1,
        "comparisons": 1,
        "blast_calls": 1,
        "sequence_bases": 30,
        "base_comparisons": 30,
        "upload_bytes": len(b">sample\nACGTTGCAAGCTTAGCGATCGATGCTAGCA\n") + len(json.dumps(custom_database()).encode()),
    }


@pytest.mark.parametrize(
    "limit,value",
    [
        ("MAX_RECORDS", 0),
        ("MAX_COMPARISONS", 0),
        ("MAX_BLAST_CALLS", 0),
        ("MAX_BASE_COMPARISONS", 27),
        ("MAX_UPLOAD_BYTES", 40),
    ],
)
@pytest.mark.parametrize("endpoint", ["/api/preflight", "/api/analyze"])
def test_preflight_and_direct_analysis_enforce_identical_limits(client, monkeypatch, limit, value, endpoint):
    monkeypatch.setattr(service, limit, value)
    monkeypatch.setattr(service, "blast_version", lambda: pytest.fail("Oversized batch launched BLAST"))
    response = custom_request(client, endpoint)
    assert response.status_code == 413, response.text
    assert response.json()["error"]["code"] in {"work_limit", "upload_too_large"}


def test_preflight_uses_shared_influenza_segment_filter(client, monkeypatch):
    monkeypatch.setattr(service, "blast_version", lambda: pytest.fail("Preflight launched BLAST"))
    database = custom_database()
    scheme = database["schemes"][0]
    scheme["organism"] = "Influenza-A"
    scheme["primers"][0].update(segment="HA", subtype_tags=["H3"])
    response = custom_request(
        client,
        "/api/preflight",
        database,
        b">01-HA|sample\nACGT\n>03-M|sample\nACGT\n",
        virus="influenza",
        flu_type="H3",
    )
    assert response.status_code == 200, response.text
    assert response.json()["workload"]["records"] == 2
    assert response.json()["workload"]["comparisons"] == 1
    assert custom_request(client, "/api/preflight", database, virus="influenza", flu_type="H1").status_code == 422


@pytest.mark.parametrize("assay_type", ["pcr", "ngs"])
def test_custom_influenza_all_segments_subtypes_and_exact_workload(client, simulated_blast, assay_type):
    database = custom_database(assay_type)
    scheme = database["schemes"][0]
    scheme["organism"] = "Influenza-A"
    original = scheme["primers"][0]
    segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS", "CUSTOM1"]
    scheme["primers"] = [
        dict(original, id=s, name=s, segment=s.lower(), subtype_tags=["h5n1", "H7N9"]) for s in segments
    ] + [
        dict(original, id="shared", name="H3_shared_PB2", segment="PB2", subtype_tags=[]),
        dict(original, id="other", name="other", segment="NA", subtype_tags=["H9N2"]),
    ]
    catalog_response = client.post("/api/database", files={"database": ("flu.json", json.dumps(database))})
    assert catalog_response.status_code == 200, catalog_response.text
    catalog = catalog_response.json()
    assert catalog["viruses"][0]["subtypes"] == ["A", "H5N1", "H7N9", "H9N2"]
    assert catalog["assays"][0]["subtype_counts"] == {"H5N1": 9, "H7N9": 9, "H9N2": 1}
    assert catalog["assays"][0]["untagged_primers"] == 1
    fasta = "".join(f">{i:02}-{s}|sample\n{original['sequence']}\n" for i, s in enumerate(segments, 1)).encode()
    selection = {"virus": "influenza", "flu_type": "H5N1", "assay_type": assay_type}
    preflight = custom_request(client, "/api/preflight", database, fasta, **selection)
    assert preflight.status_code == 200, preflight.text
    assert preflight.json()["workload"]["comparisons"] == 10
    result = custom_request(client, "/api/analyze", database, fasta, **selection)
    assert result.status_code == 200, result.text
    assert result.json()["summary"]["comparisons"] == 10
    rows = result.json()["rows"]
    assert all(row["Primer_Segment"] == row["Subject_Segment"] for row in rows)
    assert {row["Primer_Name"] for row in rows} == set(segments) | {"H3_shared_PB2"}


@pytest.mark.parametrize("tags", ["H5N1", None, [""], ["H5 N1"], ["H5/N1"], ["H" * 65]])
def test_malformed_subtype_labels_fail_validation_in_cli_and_api(client, tags):
    database = custom_database()
    database["schemes"][0]["primers"][0]["subtype_tags"] = tags
    assert engine.validate_normalized_primer_library(database).errors
    result = client.post("/api/database", files={"database": ("flu.json", json.dumps(database))})
    assert result.status_code == 422


@pytest.mark.parametrize(
    "organism,segment,tag",
    [
        ("Influenza-B", "NA", "Victoria"),
        ("Influenza-C", "HEF", "lineage1"),
        ("Influenza-D", "P3", "custom1"),
    ],
)
def test_catalog_and_preflight_support_database_defined_influenza_types(client, organism, segment, tag):
    database = custom_database()
    scheme = database["schemes"][0]
    scheme["organism"] = organism
    scheme["primers"][0].update(segment=segment, subtype_tags=[tag])
    response = client.post("/api/database", files={"database": ("flu.json", json.dumps(database))})
    assert response.status_code == 200, response.text
    selection = organism.split("-")[1] + "/" + tag.upper()
    assert response.json()["viruses"][0]["selections"][selection] == {"organism": organism, "tag": tag.upper()}
    preflight = custom_request(
        client,
        "/api/preflight",
        database,
        f">01-{segment}|sample\nACGT\n".encode(),
        virus="influenza",
        flu_type=selection,
    )
    assert preflight.status_code == 200, preflight.text
    assert preflight.json()["workload"]["comparisons"] == 1


@pytest.mark.skipif(not shutil.which("blastn"), reason="BLAST+ integration requires blastn")
def test_real_blast_h5n1_pb2_web_and_cli_agree(client, tmp_path):
    database = custom_database()
    scheme = database["schemes"][0]
    scheme["organism"] = "Influenza-A"
    scheme["primers"][0].update(segment="PB2", subtype_tags=["H5N1"])
    fasta = b">01-PB2|sample\nACGTTGCAAGCTTAGCGATCGATGCTAGCA\n>06-NA|other\nACGTTGCAAGCTTAGCGATCGATGCTAGCA\n"
    result = custom_request(client, "/api/analyze", database, fasta, virus="influenza", flu_type="H5N1")
    assert result.status_code == 200, result.text
    db_file, fasta_file, csv_file = [tmp_path / name for name in ["flu.json", "custom.fasta", "out.csv"]]
    db_file.write_text(json.dumps(database))
    fasta_file.write_bytes(fasta)
    cli = subprocess.run(
        [
            sys.executable,
            str(ROOT / "primer_checker.py"),
            "--primers",
            str(db_file),
            "--virus",
            "influenza",
            "--flu-type",
            "H5N1",
            "--assay-type",
            "pcr",
            "--fasta",
            str(fasta_file),
            "--output",
            str(csv_file),
        ],
        capture_output=True,
        text=True,
    )
    assert cli.returncode == 0, cli.stderr
    rows = list(csv.DictReader(csv_file.open()))
    assert len(rows) == len(result.json()["rows"]) == 1
    for field in [
        "Primer_Name",
        "Subject_Sequence_ID",
        "Primer_Segment",
        "Subject_Segment",
        "Mismatches",
        "Hit_Status",
    ]:
        assert rows[0][field] == str(result.json()["rows"][0][field])
