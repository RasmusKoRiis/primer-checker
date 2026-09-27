/** Original 16 × 16 pixel glyphs. Motion lives in CSS; icons are decorative. */
type Pixel = readonly [x: number, y: number, width: number, height: number];
type IconProps = { size?: number; className?: string };

function pixelIcon(
  name: string,
  main: Pixel[],
  accent: Pixel[] = [],
  detail: Pixel[] = [],
) {
  return function PixelGlyph({ size = 16, className = "" }: IconProps) {
    return (
      <svg
        width={size}
        height={size}
        viewBox="0 0 16 16"
        fill="currentColor"
        shapeRendering="crispEdges"
        aria-hidden="true"
        focusable="false"
        className={`pixel-icon pixel-${name} ${className}`}
      >
        {[main, accent, detail].map((pixels, layer) => (
          <g
            key={layer}
            className={["pixel-main", "pixel-accent", "pixel-detail"][layer]}
          >
            {pixels.map(([x, y, width, height], index) => (
              <rect key={index} x={x} y={y} width={width} height={height} />
            ))}
          </g>
        ))}
      </svg>
    );
  };
}

// Two complementary strands, with individual bases visible even at 16 px.
export const PrimerMark = pixelIcon(
  "primer",
  [
    [3, 1, 2, 4],
    [5, 5, 2, 2],
    [7, 7, 2, 2],
    [9, 9, 2, 2],
    [11, 11, 2, 4],
  ],
  [
    [11, 1, 2, 4],
    [9, 5, 2, 2],
    [7, 7, 2, 2],
    [5, 9, 2, 2],
    [3, 11, 2, 4],
  ],
  [
    [5, 2, 6, 1],
    [5, 13, 6, 1],
  ],
);

export const Database = pixelIcon(
  "database",
  [
    [3, 1, 10, 2],
    [1, 3, 2, 10],
    [13, 3, 2, 10],
    [3, 13, 10, 2],
    [3, 5, 10, 1],
    [3, 9, 10, 1],
  ],
  [
    [4, 3, 2, 2],
    [7, 6, 2, 2],
    [10, 10, 2, 2],
  ],
);

export const FileText = pixelIcon(
  "file",
  [
    [3, 1, 6, 2],
    [3, 3, 2, 10],
    [3, 13, 10, 2],
    [11, 6, 2, 7],
    [9, 1, 2, 2],
    [11, 3, 2, 2],
    [8, 3, 2, 3],
    [10, 5, 3, 1],
  ],
  [
    [6, 8, 4, 1],
    [6, 10, 4, 1],
  ],
);

export const FlaskConical = pixelIcon(
  "flask",
  [
    [5, 1, 6, 2],
    [5, 3, 2, 4],
    [9, 3, 2, 4],
    [3, 7, 2, 3],
    [11, 7, 2, 3],
    [1, 10, 2, 3],
    [13, 10, 2, 3],
    [3, 13, 10, 2],
  ],
  [
    [5, 9, 2, 2],
    [9, 10, 2, 2],
    [7, 6, 1, 1],
  ],
);

const tray: Pixel[] = [
  [2, 11, 2, 4],
  [4, 13, 8, 2],
  [12, 11, 2, 4],
];
const down: Pixel[] = [
  [7, 1, 2, 8],
  [3, 5, 2, 2],
  [5, 7, 2, 2],
  [9, 7, 2, 2],
  [11, 5, 2, 2],
  [7, 9, 2, 2],
];
const up: Pixel[] = [
  [7, 1, 2, 10],
  [5, 3, 2, 2],
  [3, 5, 2, 2],
  [9, 3, 2, 2],
  [11, 5, 2, 2],
];
export const Upload = pixelIcon("upload", tray, up);
export const Download = pixelIcon("download", tray, down);
export const ArrowDownToLine = Download;
export const ArrowDown = pixelIcon("down", [
  [7, 2, 2, 12],
  [3, 8, 2, 2],
  [5, 10, 2, 2],
  [9, 10, 2, 2],
  [11, 8, 2, 2],
]);
export const ArrowRight = pixelIcon("right", [
  [2, 7, 12, 2],
  [8, 3, 2, 2],
  [10, 5, 2, 2],
  [10, 9, 2, 2],
  [8, 11, 2, 2],
]);
export const ArrowUpRight = pixelIcon("up-right", [
  [5, 2, 9, 2],
  [12, 4, 2, 7],
  [10, 4, 2, 2],
  [8, 6, 2, 2],
  [6, 8, 2, 2],
  [4, 10, 2, 2],
  [2, 12, 2, 2],
]);
export const ChevronRight = pixelIcon("chevron", [
  [5, 2, 2, 2],
  [7, 4, 2, 2],
  [9, 6, 2, 4],
  [7, 10, 2, 2],
  [5, 12, 2, 2],
]);
export const ChevronLeft = pixelIcon("chevron", [
  [9, 2, 2, 2],
  [7, 4, 2, 2],
  [5, 6, 2, 4],
  [7, 10, 2, 2],
  [9, 12, 2, 2],
]);
export const Check = pixelIcon("check", [
  [2, 7, 2, 2],
  [4, 9, 2, 2],
  [6, 11, 2, 2],
  [8, 9, 2, 2],
  [10, 7, 2, 2],
  [12, 5, 2, 2],
  [14, 3, 2, 2],
]);
export const Plus = pixelIcon("plus", [
  [7, 2, 2, 12],
  [2, 7, 5, 2],
  [9, 7, 5, 2],
]);
export const X = pixelIcon("close", [
  [3, 3, 2, 2],
  [5, 5, 2, 2],
  [7, 7, 2, 2],
  [9, 9, 2, 2],
  [11, 11, 2, 2],
  [11, 3, 2, 2],
  [9, 5, 2, 2],
  [5, 9, 2, 2],
  [3, 11, 2, 2],
]);
export const Info = pixelIcon(
  "info",
  [
    [5, 1, 6, 2],
    [3, 3, 2, 2],
    [1, 5, 2, 6],
    [3, 11, 2, 2],
    [5, 13, 6, 2],
    [11, 11, 2, 2],
    [13, 5, 2, 6],
    [11, 3, 2, 2],
  ],
  [
    [7, 4, 2, 2],
    [7, 7, 2, 5],
  ],
);
export const ShieldCheck = pixelIcon(
  "shield",
  [
    [7, 1, 2, 1],
    [3, 2, 4, 2],
    [9, 2, 4, 2],
    [2, 4, 2, 6],
    [12, 4, 2, 6],
    [4, 10, 2, 2],
    [10, 10, 2, 2],
    [6, 12, 4, 2],
  ],
  [
    [5, 6, 2, 2],
    [7, 8, 2, 2],
    [9, 6, 2, 2],
  ],
);
export const Trash2 = pixelIcon(
  "trash",
  [
    [6, 1, 4, 2],
    [2, 3, 12, 2],
    [3, 5, 2, 8],
    [11, 5, 2, 8],
    [5, 13, 6, 2],
  ],
  [[7, 6, 2, 5]],
);
export const SlidersHorizontal = pixelIcon(
  "sliders",
  [
    [1, 3, 3, 2],
    [8, 3, 7, 2],
    [1, 11, 7, 2],
    [12, 11, 3, 2],
  ],
  [
    [4, 1, 4, 6],
    [8, 9, 4, 6],
  ],
);
export const ArrowUpDown = pixelIcon("sort", [
  [3, 2, 2, 12],
  [1, 4, 2, 2],
  [5, 4, 2, 2],
  [11, 2, 2, 12],
  [9, 10, 2, 2],
  [13, 10, 2, 2],
]);

export function PixelLoader({ size = 16, className = "" }: IconProps) {
  return (
    <svg
      width={size}
      height={size}
      viewBox="0 0 16 16"
      fill="currentColor"
      shapeRendering="crispEdges"
      aria-hidden="true"
      focusable="false"
      className={`pixel-icon pixel-loader ${className}`}
    >
      {[
        [3, 1],
        [7, 1],
        [11, 3],
        [11, 7],
        [9, 11],
        [5, 11],
        [1, 9],
        [1, 5],
      ].map(([x, y], i) => (
        <rect
          key={i}
          x={x}
          y={y}
          width="3"
          height="3"
          style={{ animationDelay: `${i * 110}ms` }}
        />
      ))}
    </svg>
  );
}
