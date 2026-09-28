import norwegian from "./norwegian.json";
export type Language = "en" | "no";
export type Params = Record<string, string | number>;
const dictionary: Record<string, string> = norwegian;
// Match messages returned by the API without translating uploaded names or data.
const patterns = Object.entries(dictionary)
  .filter(([key]) => /\{\w+\}/.test(key))
  .map(([key, translation]) => {
    const names: string[] = [];
    const source = key
      .split(/(\{\w+\})/)
      .map((part) => {
        if (/^\{\w+\}$/.test(part)) {
          names.push(part.slice(1, -1));
          return "(.+?)";
        }
        return part.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
      })
      .join("");
    return { matcher: new RegExp("^" + source + "$"), names, translation };
  });
export function translate(
  language: Language,
  key: string,
  values: Params = {},
): string {
  let template = key;
  let params = values;
  if (language === "no") {
    if (dictionary[key]) template = dictionary[key];
    else
      for (const { matcher, names, translation } of patterns) {
        const match = key.match(matcher);
        if (match) {
          template = translation;
          params = {
            ...Object.fromEntries(names.map((name, i) => [name, match[i + 1]])),
            ...values,
          };
          break;
        }
      }
  }
  return template.replace(/\{(\w+)\}/g, (match, name) =>
    String(params[name] ?? match),
  );
}
