"""
Renders the plain-text output explanations (see analysis_helper.get_pipeline_readme_content) as a
self-contained HTML page. The text is parsed instead of duplicated so that the .txt files stay the
single source of the explanations.

Recognized text layout:
    ====== / title / ======          page title ("TACTIC Pipeline — <Title>")
    Heading + line of '-'            section
    • <file name> (<note>)           file entry; consecutive '•' lines share the description below
      - What it is: ...              labeled description (continuation lines indented by 4)
            command ...              code line (indented by 6 or more)
    1. Step title + indented text    numbered step
    (...)                            note
    text + indented '•' / '1.'       paragraph with list
"""
import html
import re
from pathlib import Path
from typing import Dict, List, Optional
from urllib.parse import quote


_RULE_RE = re.compile(r"^\s*([=-])\1{3,}\s*$")
_FIELD_RE = re.compile(r"^(What (?:it is|they are)(?: used for)?):\s*(.*)$")
_NAME_NOTE_RE = re.compile(r"^(.*?\S)\s+\(((?:present|if)\b.*)\)$")
_STEP_RE = re.compile(r"^(\d+)\.\s+(.*)$")
_INNER_ITEM_RE = re.compile(r"^\s+(•|\d+\.)\s+(.*)$")
_SECTION_NUMBER_RE = re.compile(r"^(\d+)\.\s+(.*)$")
_EXPLANATION_REF_RE = re.compile(r"((?:[\w.-]+/)?)Output files explanation\.(?:txt|html)")
_FILE_NAME_RE = re.compile(r"(?<![\w/.-])([\w*.-]+\.(?:tab|fasta|tre|html|txt|csv|fastq))\b")
_FIELD_LABELS = {
    "What it is": "What it is",
    "What they are": "What they are",
    "What it is used for": "Use it for",
    "What they are used for": "Use them for",
}
_FILE_TYPES = {
    "folder": "Folder",
    "map": "Mapping",
    "table": "Table",
    "seq": "Sequences",
    "tree": "Tree",
    "chart": "Interactive",
    "log": "Log",
    "doc": "Document",
}
HTML_EXPLANATION_NAME = "Output files explanation.html"


def _is_rule(line: str, char: str) -> bool:
    match = _RULE_RE.match(line)
    return bool(match) and match.group(1) == char


def _join(lines: List[str]) -> str:
    return re.sub(r"\s{2,}", " ", " ".join(line.strip() for line in lines)).strip()


def _parse_file_entry(lines: List[str]) -> Dict:
    names, fields = [], []
    for line in lines:
        stripped = line.strip()
        indent = len(line) - len(line.lstrip())
        if line.startswith("• "):
            names.append(stripped[2:].strip())
        elif stripped.startswith("- ") and indent <= 3:
            field_match = _FIELD_RE.match(stripped[2:])
            fields.append({
                "label": _FIELD_LABELS.get(field_match.group(1), field_match.group(1)) if field_match else "",
                "text": field_match.group(2) if field_match else stripped[2:],
                "code": [],
            })
        elif fields and indent >= 6:
            fields[-1]["code"].append(stripped)
        elif fields:
            fields[-1]["text"] = _join([fields[-1]["text"], stripped])
        elif names:
            names[-1] = _join([names[-1], stripped])
    files = []
    for name in names:
        note_match = _NAME_NOTE_RE.match(name)
        if note_match:
            note = note_match.group(2)
            files.append({"name": note_match.group(1), "note": note[:1].upper() + note[1:]})
        else:
            files.append({"name": name, "note": ""})
    return {"kind": "files", "files": files, "fields": fields}


def _parse_chunk(lines: List[str]) -> Dict:
    first = lines[0]
    if first.startswith("• "):
        return _parse_file_entry(lines)
    step_match = _STEP_RE.match(first)
    if step_match and len(lines) > 1 and all(line.startswith("   ") for line in lines[1:]):
        return {"kind": "step", "number": step_match.group(1), "title": step_match.group(2), "text": _join(lines[1:])}
    if first.strip().startswith("(") and lines[-1].strip().endswith(")"):
        return {"kind": "note", "text": _join(lines)[1:-1]}
    lead, items, ordered = [], [], False
    for line in lines:
        item_match = _INNER_ITEM_RE.match(line)
        if item_match:
            ordered = item_match.group(1) != "•"
            items.append(_join([item_match.group(2)]))
        elif items and line.startswith("    "):
            items[-1] = _join([items[-1], line])
        else:
            lead.append(line)
    return {"kind": "para", "text": _join(lead), "items": items, "ordered": ordered}


def parse_readme_text(content: str) -> Dict:
    """
    Parses an output explanation text into {"title", "eyebrow", "sections": [{"heading", "blocks"}]}.
    """
    lines = content.splitlines()
    eq_rules = [ind for ind, line in enumerate(lines) if _is_rule(line, "=")]
    title, start, end = "", 0, len(lines)
    if len(eq_rules) >= 2:
        title = _join(lines[eq_rules[0] + 1:eq_rules[1]])
        start = eq_rules[1] + 1
        end = eq_rules[-1] if len(eq_rules) > 2 else len(lines)
    eyebrow, _, short_title = title.rpartition(" — ")

    sections = [{"heading": "", "chunks": []}]
    chunk: List[str] = []
    body = lines[start:end]
    ind = 0
    while ind < len(body):
        line = body[ind]
        next_line = body[ind + 1] if ind + 1 < len(body) else ""
        if line.strip() and _is_rule(next_line, "-"):
            if chunk:
                sections[-1]["chunks"].append(chunk)
                chunk = []
            sections.append({"heading": line.strip(), "chunks": []})
            ind += 2
            continue
        if line.strip():
            chunk.append(line.rstrip())
        elif chunk:
            sections[-1]["chunks"].append(chunk)
            chunk = []
        ind += 1
    if chunk:
        sections[-1]["chunks"].append(chunk)

    parsed_sections = []
    for section in sections:
        if not section["heading"] and not section["chunks"]:
            continue
        parsed_sections.append({
            "heading": section["heading"],
            "blocks": [_parse_chunk(chunk_lines) for chunk_lines in section["chunks"]],
        })
    return {"title": short_title or title, "eyebrow": eyebrow, "sections": parsed_sections}


def _inline(text: str, folder: Optional[Path] = None) -> str:
    """
    Escapes text and adds light formatting: `code`, file names, arrows and links to other explanation pages.
    If folder is given, links to explanation pages of subfolders are only added when those pages exist.
    """
    placeholders: List[str] = []

    def _keep(fragment: str) -> str:
        placeholders.append(fragment)
        return f"\x00{len(placeholders) - 1}\x00"

    def _explanation_ref(match: re.Match) -> str:
        sub_dir = match.group(1)
        if folder is not None and sub_dir and not (folder / sub_dir / HTML_EXPLANATION_NAME).is_file():
            return _keep(f'<code class="fn">{sub_dir}Output files explanation</code>')
        return _keep(f'<a href="{quote(sub_dir)}{quote(HTML_EXPLANATION_NAME)}">{sub_dir}Output files explanation</a>')

    escaped = html.escape(text, quote=False)
    escaped = _EXPLANATION_REF_RE.sub(_explanation_ref, escaped)
    escaped = re.sub(r"`([^`]+)`", lambda m: _keep(f"<code>{m.group(1)}</code>"), escaped)
    escaped = _FILE_NAME_RE.sub(lambda m: _keep(f'<code class="fn">{m.group(1)}</code>'), escaped)
    escaped = escaped.replace(" -&gt; ", " → ")
    return re.sub(r"\x00(\d+)\x00", lambda m: placeholders[int(m.group(1))], escaped)


def _file_type(name: str) -> str:
    lower = name.lower()
    if lower.endswith("/"):
        return "folder"
    if lower.startswith("map-"):
        return "map"
    if lower.endswith((".tab", ".csv")):
        return "table"
    if lower.endswith(".fasta"):
        return "seq"
    if lower.endswith(".tre"):
        return "tree"
    if lower.endswith(".html"):
        return "chart"
    if "log" in lower:
        return "log"
    return "doc"


def _find_in_folder(name: str, output_dir: Optional[Path]) -> List[Dict]:
    """
    Returns [{"label", "href"}] for files of output_dir matching the documented name. Placeholders like
    <Pipeline> and '*' are treated as wildcards. Folders link to their own explanation page.
    """
    if output_dir is None or not output_dir.is_dir():
        return []
    if name.endswith("/"):
        folder = output_dir / name.rstrip("/")
        if (folder / HTML_EXPLANATION_NAME).is_file():
            return [{"label": name, "href": quote(name) + quote(HTML_EXPLANATION_NAME)}]
        if folder.is_dir() and any(item.is_file() for item in folder.iterdir()):
            return [{"label": name, "href": quote(name)}]
        return []
    pattern = re.sub(r"<[^>]+>", "*", name)
    if any(char in pattern for char in "*?["):
        return [{"label": p.name, "href": quote(p.name)} for p in sorted(output_dir.glob(pattern)) if p.is_file()]
    return [{"label": name, "href": quote(name)}] if (output_dir / name).is_file() else []


def _render_files(block: Dict, output_dir: Optional[Path]) -> str:
    file_type = _file_type(block["files"][0]["name"]) if block["files"] else "doc"
    search_text = " ".join(
        [f["name"] for f in block["files"]] + [f["text"] for f in block["fields"]]
    ).lower()
    name_rows = []
    for file_entry in block["files"]:
        matches = _find_in_folder(file_entry["name"], output_dir)
        is_pattern = any(char in file_entry["name"] for char in "*<")
        open_links = ""
        if matches and not is_pattern:
            open_links = f'<a class="open" href="{matches[0]["href"]}">Open</a>'
        elif matches:
            open_links = '<ul class="matches">' + "".join(
                f'<li><a href="{m["href"]}">{html.escape(m["label"])}</a></li>' for m in matches
            ) + "</ul>"
        note = f'<span class="pill">{html.escape(file_entry["note"])}</span>' if file_entry["note"] else ""
        name_rows.append(
            f'<div class="name-row"><code class="fname">{html.escape(file_entry["name"])}</code>'
            f'{note}{open_links}</div>'
        )
    field_rows = []
    unlabeled = []
    for field in block["fields"]:
        code = "".join(f"<pre><code>{html.escape(code_line)}</code></pre>" for code_line in field["code"])
        if field["label"]:
            field_rows.append(f'<dt>{html.escape(field["label"])}</dt><dd>{_inline(field["text"], output_dir)}{code}</dd>')
        else:
            unlabeled.append(f"<li>{_inline(field['text'], output_dir)}{code}</li>")
    details = ""
    if field_rows:
        details += f'<dl>{"".join(field_rows)}</dl>'
    if unlabeled:
        details += f'<ul class="plain">{"".join(unlabeled)}</ul>'
    return (
        f'<article class="file" data-search="{html.escape(search_text)}">'
        f'<span class="type t-{file_type}">{_FILE_TYPES[file_type]}</span>'
        f'<div class="file-body"><div class="names">{"".join(name_rows)}</div>{details}</div>'
        f'</article>'
    )


def _render_block(block: Dict, output_dir: Optional[Path]) -> str:
    if block["kind"] == "files":
        return _render_files(block, output_dir)
    if block["kind"] == "step":
        return (
            f'<div class="step"><span class="step-no">{html.escape(block["number"])}</span>'
            f'<div><h3>{_inline(block["title"], output_dir)}</h3><p>{_inline(block["text"], output_dir)}</p></div></div>'
        )
    if block["kind"] == "note":
        return f'<p class="note">{_inline(block["text"], output_dir)}</p>'
    out = f'<p>{_inline(block["text"], output_dir)}</p>' if block["text"] else ""
    if block["items"]:
        tag = "ol" if block["ordered"] else "ul"
        out += f"<{tag}>" + "".join(f"<li>{_inline(item, output_dir)}</li>" for item in block["items"]) + f"</{tag}>"
    return out


def _slug(text: str, used: set) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-") or "section"
    candidate, counter = slug, 2
    while candidate in used:
        candidate = f"{slug}-{counter}"
        counter += 1
    used.add(candidate)
    return candidate


_CSS = """
:root{--bg:#f6f7f9;--surface:#ffffff;--surface-2:#f0f2f5;--text:#1d2330;--muted:#5b6475;--border:#e1e5eb;
--accent:#1f6feb;--accent-soft:#e7f0fe;--pill:#fff4e0;--pill-text:#8a5300;--code:#eef1f5;
--t-table:#1f6feb;--t-seq:#7c3aed;--t-tree:#15803d;--t-chart:#c2410c;--t-map:#0e7490;--t-folder:#a16207;
--t-log:#475569;--t-doc:#475569;--shadow:0 1px 2px rgba(16,24,40,.06),0 1px 3px rgba(16,24,40,.08)}
@media (prefers-color-scheme:dark){:root{--bg:#0f1218;--surface:#171b23;--surface-2:#1e2430;--text:#e6e9ef;
--muted:#9aa4b5;--border:#2a3140;--accent:#5b9dff;--accent-soft:#1a2740;--pill:#3a2a10;--pill-text:#f5c27a;
--code:#232a37;--t-table:#5b9dff;--t-seq:#a78bfa;--t-tree:#4ade80;--t-chart:#fb923c;--t-map:#22d3ee;
--t-folder:#facc15;--t-log:#94a3b8;--t-doc:#94a3b8;--shadow:none}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--text);font:16px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",
Roboto,"Helvetica Neue",Arial,sans-serif}
a{color:var(--accent)}
code,pre{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,"Liberation Mono",monospace;font-size:.88em}
code{background:var(--code);padding:.1em .35em;border-radius:4px;overflow-wrap:anywhere}
pre{background:var(--code);padding:.7em .9em;border-radius:8px;overflow-x:auto;margin:.6em 0 0}
pre code{background:none;padding:0}
header{background:var(--surface);border-bottom:1px solid var(--border)}
.wrap{max-width:1120px;margin:0 auto;padding:0 16px}
header .wrap{padding-top:28px;padding-bottom:24px}
.eyebrow{color:var(--accent);font-weight:600;font-size:.85rem;letter-spacing:.06em;text-transform:uppercase;margin:0}
h1{font-size:clamp(1.6rem,3.5vw,2.2rem);line-height:1.2;margin:.25em 0 .2em}
.folder{color:var(--muted);margin:0;font-size:.95rem}
.search{margin-top:18px;position:relative;max-width:520px}
.search input{width:100%;font:inherit;padding:.65em .9em;border:1px solid var(--border);border-radius:10px;
background:var(--surface-2);color:var(--text)}
.search input:focus{outline:2px solid var(--accent);outline-offset:1px}
.layout{display:grid;grid-template-columns:220px minmax(0,1fr);gap:32px;padding-top:28px;padding-bottom:48px}
nav.toc{position:sticky;top:20px;align-self:start;font-size:.92rem}
nav.toc p{margin:0 0 .5em;color:var(--muted);font-weight:600;font-size:.8rem;text-transform:uppercase;letter-spacing:.05em}
nav.toc ol{list-style:none;margin:0;padding:0;border-left:2px solid var(--border)}
nav.toc a{display:block;padding:.3em .8em;color:var(--muted);text-decoration:none;margin-left:-2px;
border-left:2px solid transparent}
nav.toc a:hover{color:var(--text);border-left-color:var(--accent)}
section{margin-bottom:36px;scroll-margin-top:16px}
h2{font-size:1.3rem;margin:0 0 14px;display:flex;align-items:center;gap:10px}
.sec-no{display:inline-grid;place-items:center;min-width:1.8em;height:1.8em;border-radius:50%;
background:var(--accent-soft);color:var(--accent);font-size:.85em}
section.intro{background:var(--surface);border:1px solid var(--border);border-left:4px solid var(--accent);
border-radius:12px;padding:20px 22px;box-shadow:var(--shadow)}
section.intro h2{font-size:1.15rem}
section.intro p:last-child{margin-bottom:0}
.file{display:grid;grid-template-columns:auto minmax(0,1fr);gap:14px;background:var(--surface);
border:1px solid var(--border);border-radius:12px;padding:16px 18px;margin-bottom:12px;box-shadow:var(--shadow)}
.type{font-size:.72rem;font-weight:700;letter-spacing:.04em;text-transform:uppercase;padding:.3em .6em;
border-radius:6px;align-self:start;color:var(--surface);white-space:nowrap}
.t-table{background:var(--t-table)}.t-seq{background:var(--t-seq)}.t-tree{background:var(--t-tree)}
.t-chart{background:var(--t-chart)}.t-map{background:var(--t-map)}.t-folder{background:var(--t-folder)}
.t-log{background:var(--t-log)}.t-doc{background:var(--t-doc)}
.name-row{display:flex;flex-wrap:wrap;align-items:center;gap:8px;margin-bottom:4px}
.fname{font-weight:600;font-size:.95em;background:none;padding:0}
.pill{background:var(--pill);color:var(--pill-text);font-size:.78rem;padding:.15em .6em;border-radius:999px}
.open{font-size:.82rem;font-weight:600;text-decoration:none;border:1px solid var(--accent);border-radius:999px;
padding:.05em .7em}
.open:hover{background:var(--accent-soft)}
.matches{margin:.2em 0 0;padding-left:1.1em;font-size:.88rem;width:100%}
dl{display:grid;grid-template-columns:max-content minmax(0,1fr);gap:6px 14px;margin:10px 0 0}
dt{color:var(--muted);font-weight:600;font-size:.9rem}
dd{margin:0}
ul.plain{margin:8px 0 0;padding-left:1.1em}
.note{color:var(--muted);font-style:italic;margin:0 0 12px}
.step{display:flex;gap:14px;background:var(--surface);border:1px solid var(--border);border-radius:12px;
padding:14px 18px;margin-bottom:10px;box-shadow:var(--shadow)}
.step-no{flex:none;display:grid;place-items:center;width:2em;height:2em;border-radius:50%;background:var(--accent);
color:var(--surface);font-weight:700}
.step h3{margin:0;font-size:1rem}
.step p{margin:.2em 0 0}
.empty{display:none;color:var(--muted)}
footer{border-top:1px solid var(--border);color:var(--muted);font-size:.9rem}
footer .wrap{padding-top:16px;padding-bottom:24px}
.hidden{display:none !important}
@media (max-width:820px){.layout{grid-template-columns:minmax(0,1fr)}nav.toc{display:none}
.file{grid-template-columns:minmax(0,1fr)}dl{grid-template-columns:minmax(0,1fr)}dt{margin-top:4px}}
"""

_JS = """
(function(){
  var input=document.getElementById('filter');
  if(!input){return;}
  var empty=document.getElementById('empty');
  input.addEventListener('input',function(){
    var q=input.value.trim().toLowerCase(), shown=0;
    document.querySelectorAll('article.file').forEach(function(card){
      var hit=!q||card.getAttribute('data-search').indexOf(q)!==-1;
      card.classList.toggle('hidden',!hit); if(hit){shown++;}
    });
    document.querySelectorAll('section[data-has-files]').forEach(function(sec){
      sec.classList.toggle('hidden',!!q&&!sec.querySelector('article.file:not(.hidden)'));
    });
    document.querySelectorAll('section:not([data-has-files])').forEach(function(sec){
      sec.classList.toggle('hidden',!!q);
    });
    empty.style.display=(q&&!shown)?'block':'none';
  });
  document.addEventListener('keydown',function(e){
    if(e.key==='/'&&document.activeElement!==input){e.preventDefault();input.focus();}
  });
})();
"""


def render_readme_html(content: str, output_dir: Optional[str] = None, text_file_name: str = "") -> str:
    """
    Returns a self-contained HTML page for an output explanation text. If output_dir is given, documented
    files found in it get an "Open" link.
    """
    doc = parse_readme_text(content)
    folder = Path(output_dir).resolve() if output_dir else None
    used_ids: set = set()
    toc, sections_html = [], []
    has_files = False
    for section in doc["sections"]:
        blocks_html = "".join(_render_block(block, folder) for block in section["blocks"])
        section_has_files = any(block["kind"] == "files" for block in section["blocks"])
        has_files = has_files or section_has_files
        heading = section["heading"]
        number_match = _SECTION_NUMBER_RE.match(heading)
        is_intro = bool(heading) and not number_match and not section_has_files and heading.endswith("?")
        attrs = ' data-has-files="1"' if section_has_files else ""
        if is_intro:
            attrs += ' class="intro"'
        if heading:
            sec_id = _slug(heading, used_ids)
            title_html = (
                f'<span class="sec-no">{number_match.group(1)}</span>{html.escape(number_match.group(2))}'
                if number_match else html.escape(heading)
            )
            toc.append(f'<li><a href="#{sec_id}">{html.escape(number_match.group(2) if number_match else heading)}</a></li>')
            sections_html.append(
                f'<section id="{sec_id}"{attrs}>'
                f'<h2>{title_html}</h2>{blocks_html}</section>'
            )
        else:
            sections_html.append(f'<section class="intro"{attrs}>{blocks_html}</section>')

    title = doc["title"] or "Output files explanation"
    eyebrow = doc["eyebrow"] or "TACTIC Pipeline"
    folder_line = f'<p class="folder">Folder: <code>{html.escape(folder.name)}/</code></p>' if folder else ""
    search = (
        '<div class="search"><input id="filter" type="search" '
        'placeholder="Filter files by name or description (press /)" aria-label="Filter files"></div>'
        if has_files else ""
    )
    toc_html = f'<nav class="toc"><p>On this page</p><ol>{"".join(toc)}</ol></nav>' if toc else "<div></div>"
    text_link = (
        f'Plain-text version: <a href="{quote(text_file_name)}">{html.escape(text_file_name)}</a>'
        if text_file_name else "Generated by the TACTIC pipeline."
    )
    return (
        "<!DOCTYPE html>\n"
        '<html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        f"<title>{html.escape(title)}</title><style>{_CSS}</style></head><body>"
        f'<header><div class="wrap"><p class="eyebrow">{html.escape(eyebrow)}</p>'
        f"<h1>{html.escape(title)}</h1>{folder_line}{search}</div></header>"
        f'<div class="wrap layout">{toc_html}<main>{"".join(sections_html)}'
        '<p id="empty" class="empty">No files match your filter.</p></main></div>'
        f'<footer><div class="wrap">{text_link}</div></footer>'
        f"<script>{_JS}</script></body></html>\n"
    )
