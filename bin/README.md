# BLAST deployment bundle

`scripts/package_blast.py` populates this directory during the Linux Vercel
install step. Executables, libraries, and generated manifests are gitignored.
The official archive and extracted binary both have pinned SHA-256 checksums.
Only `blastn` and resolved non-glibc dependencies are included. See
`docs/deployment.md` for Linux verification and deployment prerequisites.
