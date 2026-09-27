# Website visual style

The interface uses navy ink, cobalt actions, teal accents and warm amber on a
light neutral background. Shared colors live in `app/globals.css`. Success,
warning and error colors remain distinct from navigation and selection colors.

The icons in `app/components/pixel-icons.tsx` are original 16 × 16 glyphs made
from square SVG rectangles, animated with CSS. They replace the external icon
library. Use 16 px for controls, 20 px for section headings and 24 px for the
Primer Checker mark. The upload target uses a 32 px icon. The matching static
browser-tab icon is `app/icon.svg`.

Keep motion small and stepped. The logo plays a short entrance sequence; hover
and keyboard focus trigger brief interactions. Only the loading indicator loops
while work is pending. Under `prefers-reduced-motion: reduce`, icons and the
intro remain still. Icons are hidden from assistive technology; the surrounding
controls provide their accessible names.

Keep the opening section text-only and the working area compact. Maintain the
English/Norwegian layout, visible dummy-data notices, and readable result colors.
