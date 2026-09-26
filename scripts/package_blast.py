#!/usr/bin/env python3
"""Package the official pinned Linux BLAST binary and its non-glibc libraries.

Run during Vercel's install step, before function file collection. Runtime
downloads are intentionally unsupported. Only one verified tar member is read.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import tarfile
import tempfile
import urllib.request
from pathlib import Path

VERSION = "2.15.0"
URL = f"https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/{VERSION}/ncbi-blast-{VERSION}+-x64-linux.tar.gz"
SHA256 = "c5c0b7029069fe5b9cd5913ae606c7f2aa901db8264ca92cfaf04e4d4768c054"
BINARY_SHA256 = "6c7532867275420a52eed9599d9aff478a2d93875cd6b1294f0c5676a6b04c04"
ROOT = Path(__file__).resolve().parents[1]
# glibc and its loader must come from the execution environment, never a bundle.
SYSTEM_LIBRARIES = {"libc.so.6", "libm.so.6", "libpthread.so.0", "libdl.so.2", "libresolv.so.2", "librt.so.1"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def extract_binary(archive: Path, target: Path):
    if sha256(archive) != SHA256:
        raise RuntimeError("BLAST archive checksum mismatch; refusing to package it.")
    with tarfile.open(archive, "r:gz") as package:
        member = package.getmember(f"ncbi-blast-{VERSION}+/bin/blastn")
        if not member.isfile() or member.size > 50_000_000:
            raise RuntimeError("Unexpected BLAST archive member.")
        with package.extractfile(member) as source, target.open("wb") as destination:
            shutil.copyfileobj(source, destination)
    if sha256(target) != BINARY_SHA256:
        raise RuntimeError("BLAST binary checksum mismatch.")
    target.chmod(0o755)


def package(output: Path, archive: Path | None = None):
    if platform.system() != "Linux" or platform.machine() not in {"x86_64", "amd64"}:
        raise RuntimeError("BLAST deployment packaging requires Linux x86_64. On macOS use BLAST+ on PATH.")
    output.mkdir(parents=True, exist_ok=True)
    binary = output / "blastn"
    if not binary.exists() or sha256(binary) != BINARY_SHA256:
        with tempfile.TemporaryDirectory(prefix="blast-build-") as temp:
            source = archive or Path(temp) / "blast.tar.gz"
            if archive is None:
                print(f"Downloading official BLAST+ {VERSION} for Linux x86_64…", flush=True)
                with urllib.request.urlopen(URL, timeout=120) as response, source.open("wb") as target:
                    shutil.copyfileobj(response, target)
            extract_binary(source, binary)
    libraries = output / "lib"
    libraries.mkdir(exist_ok=True)
    environment = {**os.environ, "LD_LIBRARY_PATH": str(libraries)}
    dependencies = subprocess.run(
        ["ldd", str(binary)], check=True, capture_output=True, text=True, env=environment
    ).stdout
    if "not found" in dependencies:
        raise RuntimeError("BLAST needs libraries missing from this Linux build image:\n" + dependencies)
    for name, location in re.findall(r"^\s*(\S+) => (/\S+)", dependencies, re.MULTILINE):
        if name not in SYSTEM_LIBRARIES and Path(location).resolve() != (libraries / name).resolve():
            shutil.copy2(Path(location).resolve(), libraries / name)
    check = subprocess.run([str(binary), "-version"], check=True, capture_output=True, text=True, env=environment)
    if f"{VERSION}+" not in check.stdout:
        raise RuntimeError("Unexpected BLAST version.")
    with tempfile.TemporaryDirectory(prefix="blast-smoke-") as temp:
        fasta = Path(temp) / "smoke.fasta"
        fasta.write_text(">synthetic\nCTGCAGATTTGGATGATTTCTCC\n")
        hit = subprocess.run(
            [str(binary), "-query", str(fasta), "-subject", str(fasta), "-word_size", "4", "-outfmt", "6"],
            capture_output=True,
            text=True,
            check=True,
            timeout=10,
            env=environment,
        )
        if not hit.stdout.strip():
            raise RuntimeError("Packaged BLAST smoke analysis returned no hit.")
    files = [binary, *sorted(libraries.iterdir())]
    manifest = {
        "version": VERSION,
        "url": URL,
        "archive_sha256": SHA256,
        "files": {str(p.relative_to(output)): {"sha256": sha256(p), "bytes": p.stat().st_size} for p in files},
        "total_bytes": sum(p.stat().st_size for p in files),
    }
    if manifest["total_bytes"] > 60_000_000:
        raise RuntimeError("BLAST bundle unexpectedly exceeds 60 MB.")
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Verified BLAST bundle: {manifest['total_bytes'] / 1_000_000:.1f} MB; {len(files) - 1} shared libraries.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--archive", type=Path, help="Use a previously downloaded archive (checksum is still required)."
    )
    parser.add_argument("--output", type=Path, default=ROOT / "bin")
    args = parser.parse_args()
    package(args.output.resolve(), args.archive)
