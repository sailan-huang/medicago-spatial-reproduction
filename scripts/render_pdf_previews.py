from __future__ import annotations

import shutil
from pathlib import Path

import fitz

ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = ROOT / "results" / "figures"
OUTPUT_DIR = ROOT / "build" / "docx_assets"


def render_previews() -> list[Path]:
    if OUTPUT_DIR.exists():
        shutil.rmtree(OUTPUT_DIR)
    OUTPUT_DIR.mkdir(parents=True)
    rendered: list[Path] = []
    for pdf_path in sorted(FIGURE_DIR.glob("*.pdf")):
        with fitz.open(pdf_path) as document:
            if len(document) != 1:
                raise ValueError(f"Expected one-page figure: {pdf_path}")
            pixmap = document[0].get_pixmap(matrix=fitz.Matrix(2.4, 2.4), alpha=False)
            output = OUTPUT_DIR / f"{pdf_path.stem}.png"
            pixmap.save(output)
            rendered.append(output)
    if not rendered:
        raise FileNotFoundError(f"No PDF figures found in {FIGURE_DIR}")
    return rendered


if __name__ == "__main__":
    for path in render_previews():
        print(path)
