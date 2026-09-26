import type { Metadata } from "next";
import "./globals.css";
import { LanguageProvider } from "./components/language";

export const metadata: Metadata = {
  title: "Primer Checker · Sequence compatibility",
  description:
    "Evaluate PCR and sequencing primer compatibility against viral consensus sequences.",
};
export default function RootLayout({
  children,
}: Readonly<{ children: React.ReactNode }>) {
  return (
    <html lang="en">
      <body>
        <LanguageProvider>{children}</LanguageProvider>
      </body>
    </html>
  );
}
