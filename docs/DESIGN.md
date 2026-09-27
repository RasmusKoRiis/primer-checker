# Website visual style

The interface takes inspiration from the **Sunny Beach Day** palette on
[Coolors Trending Palettes](https://coolors.co/palettes/trending): slate
`#264653`, teal `#2A9D8F`, sand `#E9C46A`, apricot `#F4A261` and coral `#E76F51`.
The website uses a restrained selection of these colors on warm ivory surfaces.
Shared color roles live in `app/globals.css`.

- Slate ink and warm paper keep the scientific workspace readable.
- Darker teal `#21776D` marks actions, links and selected controls.
- Darker coral `#B94F39` accents the introduction and small pixel icons.
- Bright teal and coral are decorative accents; sand highlights the logo and
  mismatch summaries. The browser-tab icon uses the same slate and sand.
- Green success, amber warnings and red errors retain their separate meanings.

Use the darker variants for text instead of the original bright swatches. Teal
on paper has a contrast ratio of 5.30:1, coral on the page background 4.54:1,
and muted text on the page background 4.78:1. Preserve at least 4.5:1 for normal
text when adjusting this palette. Do not use color alone to explain a result.

The icons in `app/components/pixel-icons.tsx` are original 16 × 16 glyphs made
from square SVG rectangles, animated with CSS. They replace the external icon
library. Use 16 px for controls, 20 px for section headings and 24 px for the
Primer Checker mark. The upload target uses a 32 px icon. The matching static
browser-tab icon is `app/icon.svg`.

Keep motion small and stepped. The logo plays a short entrance sequence; hover
and keyboard focus trigger brief interactions. Loading indicators loop only
while work is pending. Under `prefers-reduced-motion: reduce`, icons and the
intro remain still. Icons are hidden from assistive technology; the surrounding
controls provide their accessible names.

The Analyze sequences button contains an original 40 × 36 pixel thermal cycler
from `app/components/pcr-machine.tsx`, displayed at 60 × 54 px. Its lid closes
in six steps when a request begins, then the display and indicator lamp cycle
until that request finishes or fails. Idle shows the open lid and tube block.
This is a decorative activity indicator, with no invented temperatures, cycle
counts or progress estimates. It adds no media download, timer or delay to the
analysis. Reduced-motion users get a static closed lid and lit display while
the translated button text and live status explain that analysis is running.

Keep the opening section text-only and the working area compact. Maintain the
English/Norwegian layout, visible dummy-data notices, and readable result colors.
