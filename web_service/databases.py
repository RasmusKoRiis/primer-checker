"""Bounded, self-contained upload formats; never resolve uploaded asset paths."""

import json
from typing import Annotated, Literal

from pydantic import BaseModel, Field, StringConstraints, ValidationError

import primer_analysis as engine

MAX_DATABASE_BYTES = 250_000
MAX_PRIMERS = 500
MAX_PRIMER_LENGTH = 200
Text = Annotated[str, StringConstraints(strip_whitespace=True, min_length=1, max_length=200)]
OptionalText = Annotated[str, StringConstraints(strip_whitespace=True, max_length=200)]


class Primer(BaseModel):
    id: Text
    name: Text
    sequence: Annotated[str, StringConstraints(strip_whitespace=True, min_length=1, max_length=MAX_PRIMER_LENGTH)]
    role: Text
    segment: OptionalText
    pool: OptionalText = ""
    strand: Literal["", "+", "-"] = ""
    subtype_tags: list[
        Annotated[str, StringConstraints(strip_whitespace=True, pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")]
    ] = Field(default_factory=list, max_length=32)


class Scheme(BaseModel):
    scheme_id: Text
    display_name: Text
    organism: Text
    version: Text
    assay_type: Literal["pcr", "ngs"] = "pcr"
    primers: list[Primer] = Field(min_length=1, max_length=MAX_PRIMERS)


class Database(BaseModel):
    schema_version: Literal["1.0"]
    database_version: Text
    schemes: list[Scheme] = Field(min_length=1, max_length=MAX_PRIMERS)


def unique_keys(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("Database JSON contains duplicate keys.")
        result[key] = value
    return result


def load_uploaded(data: bytes):
    """Validate types first, then use the CLI's validation and record conversion."""
    try:
        raw = json.loads(data.decode("utf-8-sig"), object_pairs_hook=unique_keys)
    except (UnicodeDecodeError, json.JSONDecodeError, RecursionError):
        raise ValueError("The database must be a valid UTF-8 JSON file.") from None
    if not isinstance(raw, dict) or not raw:
        raise ValueError("The database must be a nonempty JSON object.")

    def check_strings(value):
        if isinstance(value, str) and any(ord(c) < 32 or ord(c) == 127 for c in value):
            raise ValueError("Database fields must not contain control characters or line breaks.")
        if isinstance(value, dict):
            for key, item in value.items():
                check_strings(key)
                check_strings(item)
        elif isinstance(value, list):
            for item in value:
                check_strings(item)

    try:
        check_strings(raw)
    except RecursionError:
        raise ValueError("Database JSON is nested too deeply.") from None
    if "panels" in raw or "viruses" in raw:
        raise ValueError(
            "Upload self-contained schemes with primer sequences, or a legacy primer dictionary. Asset-backed panels are supported through the CLI."
        )
    if "schemes" not in raw:
        # Legacy uploads are also bounded before passing through the shared loader.
        for organism, primers in raw.items():
            if not isinstance(primers, dict) or not primers or len(organism) > 200:
                raise ValueError("Legacy databases must map organism names to primer name/sequence objects.")
            for name, sequence in primers.items():
                if len(name) > 200 or not isinstance(sequence, str) or len(sequence.strip()) > MAX_PRIMER_LENGTH:
                    raise ValueError(f"Primer names and sequences must have at most {MAX_PRIMER_LENGTH} characters.")
        validation = engine.validate_legacy_primer_library(raw)
        if validation.errors:
            raise ValueError(validation.errors[0])
        records = engine.legacy_library_to_records(raw)
    else:
        try:
            normalized = Database.model_validate(raw).model_dump(exclude_unset=True)
        except ValidationError as exc:
            first = exc.errors(include_input=False, include_url=False)[0]
            location = ".".join(str(part) for part in first["loc"])
            raise ValueError(f"Database {location}: {first['msg']}.") from None
        validation = engine.validate_normalized_primer_library(normalized)
        if validation.errors:
            raise ValueError(validation.errors[0])
        records = engine.normalized_library_to_records(normalized)
    for group in records.values():
        names = [(p.scheme_id, p.name) for p in group]
        if len(names) != len(set(names)):
            raise ValueError("Primer names must be distinct within each scheme.")
    if sum(map(len, records.values())) > MAX_PRIMERS:
        raise ValueError(f"A custom database may contain at most {MAX_PRIMERS} primers.")
    if any(name.casefold() == "influenza" for name in records):
        raise ValueError(
            "Specify the influenza type in organism, for example Influenza-A or Influenza-B; put subtype labels in subtype_tags."
        )
    return records, validation
