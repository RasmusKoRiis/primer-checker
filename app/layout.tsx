import type { Metadata } from "next";
import "./globals.css";

export const metadata: Metadata = {
  title: "Primer Checker · Sequence compatibility",
  description: "Evaluate diagnostic and sequencing primer compatibility against viral consensus sequences.",
};
export default function RootLayout({ children }: Readonly<{ children: React.ReactNode }>) {
  return <html lang="en"><body>{children}</body></html>;
}
