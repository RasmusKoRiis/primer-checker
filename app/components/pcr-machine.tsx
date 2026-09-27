/** Original pixel thermal cycler. The request state drives all motion in CSS. */
export default function PcrMachine({ running }: { running: boolean }) {
  return (
    <svg
      width="60"
      height="54"
      viewBox="0 0 40 36"
      shapeRendering="crispEdges"
      aria-hidden="true"
      focusable="false"
      className="pcr-machine"
      data-running={running}
    >
      <path className="pcr-shadow" d="M4 33h32v2H4z" />
      <path className="pcr-outline" d="M8 30h4v4H8zM28 30h4v4h-4z" />
      <path className="pcr-hinge" d="M6 10h2v12H6zM32 10h2v12h-2z" />

      {/* Exposed wells and capped tubes, covered by the lid during a run. */}
      <path className="pcr-outline" d="M6 18h28v2h2v4H4v-4h2z" />
      <path className="pcr-block" d="M6 20h28v3H6z" />
      {[9, 15, 21, 27].map((x) => (
        <g key={x}>
          <rect className="pcr-outline" x={x} y="16" width="4" height="5" />
          <rect className="pcr-tube" x={x + 1} y="17" width="2" height="3" />
          <rect className="pcr-cap" x={x - 1} y="15" width="6" height="2" />
        </g>
      ))}

      <path className="pcr-outline" d="M4 22h32v2h2v8H2v-8h2z" />
      <path className="pcr-shell" d="M4 24h32v6H4z" />
      <path className="pcr-shade" d="M30 24h6v6h-6zM4 30h32v1H4z" />
      <path className="pcr-outline" d="M6 25h15v4H6z" />
      <path className="pcr-screen" d="M7 26h13v2H7z" />
      <g className="pcr-display">
        {[8, 12, 16].map((x, i) => (
          <rect
            key={x}
            x={x}
            y="26"
            width="2"
            height="2"
            style={{ animationDelay: `${i * 400}ms` }}
          />
        ))}
      </g>
      <path className="pcr-lamp" d="M25 26h2v2h-2z" />
      <path className="pcr-vent" d="M30 25h4v1h-4zM30 27h4v1h-4z" />

      <g className="pcr-lid">
        <path className="pcr-outline" d="M15 2h10v3h9v2h2v6H4V7h2V5h9z" />
        <path className="pcr-shell" d="M6 7h28v4H6z" />
        <path className="pcr-lid-band" d="M6 7h28v2H6z" />
        <path className="pcr-shade" d="M30 9h4v2h-4z" />
        <path className="pcr-handle" d="M17 3h6v2h-6z" />
        <path className="pcr-outline" d="M16 10h8v1h-8z" />
      </g>
    </svg>
  );
}
