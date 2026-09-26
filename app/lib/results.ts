export interface Catalog {
  application_version: string;
  database: {
    version: string;
    sha256: string;
    source?: string;
    filename?: string;
  };
  warnings?: string[];
  viruses: { id: string; name: string; subtypes: string[] }[];
  assays: {
    id: string;
    name: string;
    organism: string;
    type: string;
    primers: number;
  }[];
  limits: {
    upload_bytes: number;
    files: number;
    records: number;
    comparisons: number;
    blast_calls?: number;
    base_comparisons?: number;
    database_bytes?: number;
    database_primers?: number;
    primer_length?: number;
  };
}
export interface ResultRow {
  [key: string]: string | number;
  Fasta_File: string;
  Virus_Type: string;
  Assay_Type: string;
  Assay_ID: string;
  Assay_Name: string;
  Primer_Name: string;
  Primer_Sequence: string;
  Primer_Segment: string;
  Primer_Pool: string;
  Subject_Sequence_ID: string;
  Subject_Segment: string;
  Hit_Status: string;
  Percent_Identity: number | string;
  Mismatches: number | string;
  Mismatch_Positions: string;
  Mismatch_Details: string;
  Query_Alignment: string;
  Subject_Alignment: string;
  Subject_Start: number | string;
  Subject_End: number | string;
}
export interface Analysis {
  rows: ResultRow[];
  summary: {
    files: number;
    sequence_records: number;
    samples: number;
    primers: number;
    comparisons: number;
    hits: number;
    mismatch_comparisons: number;
    samples_affected: number;
  };
  manifest: {
    analysis_utc: string;
    application_version: string;
    git_commit: string;
    database: Catalog["database"];
    selection: {
      virus: string;
      flu_type: string | null;
      assay_type: string;
      assay_id: string | null;
    };
    blast: { version: string };
    files: { filename: string; records: number; sha256: string }[];
  };
  warnings: string[];
  downloads: { csv: string; html: string; manifest: string };
}
export interface Filters {
  primer: string;
  sample: string;
  segment: string;
  assay: string;
  status: string;
  mismatches: string;
}
export const emptyFilters: Filters = {
  primer: "",
  sample: "",
  segment: "",
  assay: "",
  status: "",
  mismatches: "",
};
export function filterRows(rows: ResultRow[], f: Filters) {
  return rows.filter(
    (r) =>
      r.Primer_Name.toLowerCase().includes(f.primer.toLowerCase()) &&
      r.Subject_Sequence_ID.toLowerCase().includes(f.sample.toLowerCase()) &&
      (!f.segment || (r.Primer_Segment || r.Subject_Segment) === f.segment) &&
      (!f.assay || r.Assay_ID === f.assay) &&
      (!f.status || r.Hit_Status === f.status) &&
      (!f.mismatches ||
        (r.Hit_Status === "hit" &&
          (f.mismatches === "any"
            ? Number(r.Mismatches) > 0
            : Number(r.Mismatches) === Number(f.mismatches)))),
  );
}
export function aggregatePrimers(rows: ResultRow[]) {
  const groups = new Map<
    string,
    {
      key: string;
      primer: string;
      assay: string;
      segment: string;
      pool: string;
      tested: number;
      perfect: number;
      affected: number;
      noHit: number;
      maximum: number;
      percent: number;
    }
  >();
  for (const r of rows) {
    const key = JSON.stringify([r.Assay_ID, r.Primer_Name, r.Primer_Sequence]);
    const item = groups.get(key) || {
      key,
      primer: r.Primer_Name,
      assay: r.Assay_Name || r.Assay_ID,
      segment: r.Primer_Segment,
      pool: r.Primer_Pool,
      tested: 0,
      perfect: 0,
      affected: 0,
      noHit: 0,
      maximum: 0,
      percent: 0,
    };
    item.tested++;
    if (r.Hit_Status !== "hit") item.noHit++;
    else {
      const n = Number(r.Mismatches);
      if (n === 0) item.perfect++;
      else item.affected++;
      item.maximum = Math.max(item.maximum, n);
    }
    item.percent = (100 * item.affected) / item.tested;
    groups.set(key, item);
  }
  return [...groups.values()];
}
export function validateFiles(
  files: { name: string; size: number }[],
  limits: Catalog["limits"],
  database: { name: string; size: number } | null = null,
): string | null {
  if (files.length > limits.files)
    return `Choose at most ${limits.files} FASTA files.`;
  if (files.some((f) => !/\.(fasta|fa|fas|fna)$/i.test(f.name)))
    return "Choose FASTA files (.fasta, .fa, .fas, or .fna). FASTQ and ZIP are not supported.";
  if (files.some((f) => f.size === 0))
    return "One of the FASTA files is empty.";
  if (
    files.reduce((n, f) => n + f.size, database?.size || 0) >
    limits.upload_bytes
  )
    return "Combined uploads exceed 3 MB. Split the analysis into smaller batches.";
  if (new Set(files.map((f) => f.name)).size !== files.length)
    return "Choose files with distinct filenames.";
  return null;
}
export function download(text: string, name: string, type: string) {
  const url = URL.createObjectURL(new Blob([text], { type }));
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = name;
  anchor.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
