"use client";
import {
  createContext,
  useContext,
  useEffect,
  useSyncExternalStore,
} from "react";
import { translate, type Language, type Params } from "../lib/translate";

const storageKey = "primer-checker-language";
let fallback: Language = "en";
function getLanguage(): Language {
  try {
    const saved = localStorage.getItem(storageKey);
    return saved === "no" || saved === "en" ? saved : fallback;
  } catch {
    return fallback;
  }
}
function subscribe(callback: () => void) {
  window.addEventListener("storage", callback);
  window.addEventListener("primer-language", callback);
  return () => {
    window.removeEventListener("storage", callback);
    window.removeEventListener("primer-language", callback);
  };
}
function setLanguage(language: Language) {
  fallback = language;
  try {
    localStorage.setItem(storageKey, language);
  } catch {
    /* Keep the choice for this session. */
  }
  window.dispatchEvent(new Event("primer-language"));
}
const LanguageContext = createContext({
  language: "en" as Language,
  locale: "en-GB",
  t: (key: string, values?: Params) => translate("en", key, values),
});
export function LanguageProvider({ children }: { children: React.ReactNode }) {
  const language = useSyncExternalStore(
    subscribe,
    getLanguage,
    () => "en" as Language,
  );
  useEffect(() => {
    document.documentElement.lang = language === "no" ? "nb" : "en";
    document.title = translate(
      language,
      "Primer Checker · Sequence compatibility",
    );
  }, [language]);
  return (
    <LanguageContext.Provider
      value={{
        language,
        locale: language === "no" ? "nb-NO" : "en-GB",
        t: (key, values) => translate(language, key, values),
      }}
    >
      {children}
    </LanguageContext.Provider>
  );
}
export function useLanguage() {
  return useContext(LanguageContext);
}
export function LanguageSwitch() {
  const { language, t } = useLanguage();
  return (
    <div
      className="language-switch"
      role="group"
      aria-label={t("Site language")}
    >
      <button
        type="button"
        lang="en"
        aria-pressed={language === "en"}
        onClick={() => setLanguage("en")}
      >
        English
      </button>
      <button
        type="button"
        lang="nb"
        aria-pressed={language === "no"}
        onClick={() => setLanguage("no")}
      >
        Norsk
      </button>
    </div>
  );
}
