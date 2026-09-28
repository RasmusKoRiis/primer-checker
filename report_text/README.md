# Editing report text

Edit `english.json` for English wording and `norwegian.json` for Norwegian
wording. Change only the text values on the right-hand side of each key.

Keep the JSON syntax, key names, and placeholders such as `{sample}`,
`{percent}`, `{count}`, `{visible}`, and `{total}` unchanged. Both language
files must contain the same keys.

Regenerate the HTML report after saving your changes. Existing generated HTML
files do not update automatically.

The web application's Documentation tab also imports these files directly.
Rebuild/redeploy the website after editing them. Report-specific sections are
labelled separately from the website table guide, since the two interfaces
have different filters and summary behavior. Keep calculation wording aligned
with `primer_analysis.py` and review-level rules with `primer_report.py`.
