export interface PrimerDraft {
  key: number;
  name: string;
  sequence: string;
  role: string;
  segment: string;
  pool: string;
  subtype: string;
}
export interface DatabaseDraft {
  name: string;
  organism: string;
  version: string;
  assayType: string;
  primers: PrimerDraft[];
}
export function makeDatabase(draft: DatabaseDraft): string {
  if (![draft.name, draft.organism, draft.version].every((v) => v.trim()))
    throw new Error("Enter a database name, organism, and version.");
  const names = draft.primers.map((p) => p.name.trim());
  if (new Set(names).size !== names.length)
    throw new Error("Give each primer a distinct name.");
  const primers = draft.primers.map((p, index) => {
    const sequence = p.sequence.replace(/\s/g, "").toUpperCase();
    if (!p.name.trim() || !sequence)
      throw new Error(`Enter a name and sequence for primer ${index + 1}.`);
    if (!/^[ACGTRYSWKMBDHVN]+$/.test(sequence))
      throw new Error(
        `Primer ${index + 1}: use DNA IUPAC bases only (no gaps or U).`,
      );
    return {
      id: p.name.trim(),
      name: p.name.trim(),
      sequence,
      role: p.role,
      segment: p.segment.trim().toUpperCase(),
      pool: p.pool.trim(),
      subtype_tags: [
        ...new Set(
          p.subtype
            .split(",")
            .map((tag) => tag.trim().toUpperCase())
            .filter(Boolean),
        ),
      ],
    };
  });
  return JSON.stringify(
    {
      schema_version: "1.0",
      database_version: draft.version.trim(),
      schemes: [
        {
          scheme_id:
            "custom-" +
            (draft.name
              .trim()
              .toLowerCase()
              .replace(/[^a-z0-9]+/g, "-")
              .replace(/^-|-$/g, "")
              .slice(0, 100) || "primers"),
          display_name: draft.name.trim(),
          organism: draft.organism.trim(),
          version: draft.version.trim(),
          assay_type: draft.assayType,
          primers,
        },
      ],
    },
    null,
    2,
  );
}
