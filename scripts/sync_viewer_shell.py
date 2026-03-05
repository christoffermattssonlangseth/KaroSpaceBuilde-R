#!/usr/bin/env python3
from __future__ import annotations

import argparse
import ast
from pathlib import Path
import warnings


TOKENS = {
    "title": "__KAROSPACE_TITLE__",
    "min_panel_size": "__KAROSPACE_MIN_PANEL_SIZE__",
    "max_panel_size": "__KAROSPACE_MAX_PANEL_SIZE__",
    "spot_size": "__KAROSPACE_SPOT_SIZE__",
    "data_json": "__KAROSPACE_DATA_JSON__",
    "palette_json": "__KAROSPACE_PALETTE_JSON__",
    "metadata_labels_json": "__KAROSPACE_METADATA_LABELS_JSON__",
    "outline_by_json": "__KAROSPACE_OUTLINE_BY_JSON__",
    "viewer_info_html_json": "__KAROSPACE_VIEWER_INFO_HTML_JSON__",
    "viewer_info_html": "__KAROSPACE_VIEWER_INFO_HTML__",
    "theme_icon": "__KAROSPACE_THEME_ICON__",
    "initial_theme": "__KAROSPACE_INITIAL_THEME__",
    "favicon_link": "__KAROSPACE_FAVICON_LINK__",
    "footer_logo": "__KAROSPACE_FOOTER_LOGO__",
    "background": "__KAROSPACE_BACKGROUND__",
    "text_color": "__KAROSPACE_TEXT_COLOR__",
    "header_bg": "__KAROSPACE_HEADER_BG__",
    "panel_bg": "__KAROSPACE_PANEL_BG__",
    "border_color": "__KAROSPACE_BORDER_COLOR__",
    "input_bg": "__KAROSPACE_INPUT_BG__",
    "muted_color": "__KAROSPACE_MUTED_COLOR__",
    "hover_bg": "__KAROSPACE_HOVER_BG__",
    "graph_color": "__KAROSPACE_GRAPH_COLOR__",
}


def load_html_template(exporter_path: Path) -> str:
    warnings.filterwarnings("ignore", category=SyntaxWarning)
    module = ast.parse(exporter_path.read_text(encoding="utf-8"))
    for node in module.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "HTML_TEMPLATE":
                    return ast.literal_eval(node.value)
    raise RuntimeError(f"Could not find HTML_TEMPLATE in {exporter_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Vendor the current KaroSpace viewer shell into this repo.")
    parser.add_argument(
        "--reference",
        default="/Users/chrislangseth/work/karolinska_institutet/projects/KaroSpace",
        help="Path to the KaroSpace reference repository",
    )
    parser.add_argument(
        "--output",
        default="inst/viewer/karospace_viewer_shell.html",
        help="Output HTML shell path",
    )
    args = parser.parse_args()

    reference_root = Path(args.reference).resolve()
    exporter_path = reference_root / "karospace" / "exporter.py"
    template = load_html_template(exporter_path)

    footer_logo = (
        '<div class="footer-logo">'
        '<span>KaroSpace</span>'
        '<span class="footer-link">Standalone R Export</span>'
        "</div>"
    )

    format_values = dict(TOKENS)
    format_values["footer_logo"] = footer_logo

    html = template.format(**format_values)

    output_path = Path(args.output).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")


if __name__ == "__main__":
    main()
