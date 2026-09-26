"""Vercel ASGI entrypoint; run locally with uvicorn api.index:app."""

import json

from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse, Response
from starlette.concurrency import run_in_threadpool
from starlette.datastructures import UploadFile
from starlette.exceptions import HTTPException

import primer_analysis as engine
from web_service import service

app = FastAPI(
    title="Primer Checker API",
    version=service.APP_VERSION,
    docs_url="/api/docs",
    openapi_url="/api/openapi.json",
    redoc_url=None,
)


def error(message: str, code: str, status: int):
    return JSONResponse(
        {"error": {"code": code, "message": message}},
        status_code=status,
        headers={"Cache-Control": "no-store", "X-Content-Type-Options": "nosniff"},
    )


class BodyLimit:
    """Bound the actual body, including chunked uploads, before multipart parsing."""

    def __init__(self, app):
        self.app = app

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http":
            return await self.app(scope, receive, send)
        if scope["method"] == "POST":
            chunks, size = [], 0
            while True:
                event = await receive()
                if event["type"] == "http.disconnect":
                    return
                size += len(event.get("body", b""))
                if size > service.MAX_REQUEST_BYTES:
                    return await error("Upload too large. Combined files must be under 3 MB.", "upload_too_large", 413)(
                        scope, receive, send
                    )
                chunks.append(event)
                if not event.get("more_body", False):
                    break
            iterator = iter(chunks)
            original_receive = receive

            async def bounded_receive():
                return next(iterator, None) or await original_receive()

            receive = bounded_receive
        await self.app(scope, receive, send)


app.add_middleware(BodyLimit)


@app.exception_handler(service.WebError)
async def web_error(_request, exc):
    return error(exc.message, exc.code, exc.status)


@app.exception_handler(engine.BlastError)
async def blast_error(_request, exc):
    if isinstance(exc, engine.AnalysisTimeoutError):
        return error(str(exc), "analysis_timeout", 504)
    if isinstance(exc, engine.BlastUnavailableError):
        return error("BLAST is unavailable on this deployment. Contact the maintainer.", "blast_unavailable", 503)
    return error(
        "BLAST failed to analyze these files. Check the sequences or try a smaller batch.", "blast_failure", 502
    )


@app.exception_handler(Exception)
async def unexpected_error(_request, _exc):
    return error("Analysis could not be completed. Try fewer files or contact the maintainer.", "analysis_error", 500)


@app.get("/api/catalog")
def get_catalog():
    return JSONResponse(service.catalog(), headers={"Cache-Control": "no-store"})


@app.get("/api/health")
def health():
    _, version = service.blast_version()
    service.load_database()
    return {"status": "ok", "application_version": service.APP_VERSION, "blast_version": version}


@app.post("/api/analyze")
async def analyze(request: Request):
    if not request.headers.get("content-type", "").startswith("multipart/form-data"):
        raise service.WebError("Send FASTA files as a multipart form upload.")
    try:
        async with request.form(
            max_files=service.MAX_FILES + 1, max_fields=4, max_part_size=service.MAX_UPLOAD_BYTES
        ) as form:
            allowed = {"files", "metadata", "virus", "flu_type", "assay_type", "assay_id"}
            if set(form) - allowed:
                raise service.WebError("The upload contains unsupported form fields.")
            for key in allowed - {"files"}:
                if len(form.getlist(key)) > 1:
                    raise service.WebError("Each analysis option may be supplied only once.")
            uploads = form.getlist("files")
            if not uploads or any(not isinstance(f, UploadFile) for f in uploads):
                raise service.WebError("Choose at least one FASTA file.", "missing_file")
            metadata = form.get("metadata")
            if metadata is not None and not isinstance(metadata, UploadFile):
                raise service.WebError("Metadata must be a CSV file.", "invalid_metadata")
            options = {}
            for name in ("virus", "flu_type", "assay_type", "assay_id"):
                value = form.get(name, "pcr" if name == "assay_type" else "")
                if not isinstance(value, str) or len(value) > 200:
                    raise service.WebError("Invalid analysis selection.")
                options[name] = value or None
            files = [(f.filename or "", await f.read(service.MAX_UPLOAD_BYTES + 1)) for f in uploads]
            meta = (metadata.filename or "", await metadata.read(service.MAX_UPLOAD_BYTES + 1)) if metadata else None
            result = await run_in_threadpool(service.analyze, files, meta, **options)
    except HTTPException:
        raise service.WebError("Malformed upload or too many files/form fields.", "invalid_upload") from None
    encoded = json.dumps(result, ensure_ascii=True, allow_nan=False, separators=(",", ":"))
    if len(encoded.encode()) > service.MAX_RESPONSE_BYTES:
        raise service.WebError(
            "The result exceeds the web response limit. Use fewer sequences or select one assay; larger batches can run through the CLI.",
            "result_too_large",
            413,
        )
    return Response(
        encoded,
        media_type="application/json",
        headers={"Cache-Control": "no-store", "X-Content-Type-Options": "nosniff"},
    )
