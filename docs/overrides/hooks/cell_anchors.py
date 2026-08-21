"""
MkDocs hook: add stable IDs to every Jupyter cell on notebook pages.

- Each cell (code or markdown) gets `id="cell-N"`, numbered in document order.
- Code cell output blocks additionally get `id="cell-N-output"`.
- A small anchor link is injected on hover for easy link-grabbing.
"""

import re


# Combined regex: outer cell wrappers (jp-CodeCell / jp-MarkdownCell)
# AND output wrappers (jp-Cell-outputWrapper). We dispatch in code based on
# which class is present, so a single left-to-right scan keeps order correct.
_TARGET_RE = re.compile(
    r'<div\b([^>]*?\bclass="[^"]*\b'
    r'(?:jp-CodeCell|jp-MarkdownCell|jp-Cell-outputWrapper)\b'
    r'[^"]*"[^>]*)>',
    re.IGNORECASE,
)


def _inject_id(attrs, anchor_id):
    """Return a new opening <div ...> with id injected, or None if one exists."""
    if re.search(r'\bid\s*=\s*"', attrs):
        return None
    return f'<div id="{anchor_id}"{attrs}>'


_ICON_SVG = (
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 -960 960 960" '
    'width="1em" height="1em" fill="currentColor" aria-hidden="true">'
    '<path d="M318-120q-82 0-140-58t-58-140q0-40 15-76t43-64l134-133 56 56-134 134'
    'q-17 17-25.5 38.5T200-318q0 49 34.5 83.5T318-200q23 0 45-8.5t39-25.5l133-134 57 57'
    '-134 133q-28 28-64 43t-76 15Zm79-220-57-57 223-223 57 57-223 223Zm251-28-56-57 134-133'
    'q17-17 25-38t8-44q0-50-34-85t-84-35q-23 0-44.5 8.5T558-726L425-592l-57-56 134-134'
    'q28-28 64-43t76-15q82 0 139.5 58T839-641q0 39-14.5 75T782-502L648-368Z"/>'
    '</svg>'
)


def _headerlink(anchor_id, label):
    return (
        f'<a class="cell-headerlink" href="#{anchor_id}" '
        f'title="Permanent link to {label}">{_ICON_SVG}</a>'
    )


def on_page_content(html, page, config, files, **kwargs):
    # Only act on notebook-derived pages. mkdocs-jupyter sets nb_url on these
    # when include_source: true.
    if not getattr(page, "nb_url", None):
        return html

    state = {"n": 0, "current": None}

    def _replace(match):
        attrs = match.group(1)
        class_match = re.search(r'class="([^"]*)"', attrs)
        classes = class_match.group(1).split() if class_match else []

        if "jp-CodeCell" in classes or "jp-MarkdownCell" in classes:
            state["n"] += 1
            n = state["n"]
            state["current"] = n
            replacement = _inject_id(attrs, f"cell-{n}")
            if replacement is None:
                return match.group(0)
            return replacement + _headerlink(f"cell-{n}", f"cell {n}")

        if "jp-Cell-outputWrapper" in classes and state["current"] is not None:
            n = state["current"]
            replacement = _inject_id(attrs, f"cell-{n}-output")
            if replacement is None:
                return match.group(0)
            return replacement + _headerlink(f"cell-{n}-output", f"cell {n} output")

        return match.group(0)

    return _TARGET_RE.sub(_replace, html)
