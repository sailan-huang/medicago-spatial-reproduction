from __future__ import annotations

import re
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree

ROOT = Path(__file__).resolve().parents[1]
DOCX = ROOT / "docs" / "medicago_spatial_reproduction_report.docx"
WORD_NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
DRAWING_NS = {
    "w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main",
    "a": "http://schemas.openxmlformats.org/drawingml/2006/main",
}


def validate_docx() -> None:
    if not DOCX.exists() or DOCX.stat().st_size < 200_000:
        raise AssertionError(f"DOCX missing or unexpectedly small: {DOCX}")
    with zipfile.ZipFile(DOCX) as archive:
        names = set(archive.namelist())
        required = {
            "[Content_Types].xml",
            "word/document.xml",
            "word/styles.xml",
            "word/_rels/document.xml.rels",
        }
        missing = required.difference(names)
        if missing:
            raise AssertionError(f"DOCX is missing package parts: {sorted(missing)}")
        xml = archive.read("word/document.xml")
        root = ElementTree.fromstring(xml)
        text_nodes = root.findall(".//w:t", WORD_NS)
        text = "".join(node.text or "" for node in text_nodes)
        if "苜蓿根瘤空间转录组" not in text or "结论" not in text:
            raise AssertionError("Expected Chinese report text is missing.")
        for marker in ("�", "锟斤拷", "Ã", "â€"):
            if marker in text:
                raise AssertionError(f"Possible mojibake marker found: {marker!r}")
        images = root.findall(".//a:blip", DRAWING_NS)
        if len(images) != 9:
            raise AssertionError(f"Expected 9 embedded figures, found {len(images)}.")
        if not re.search(r"<w:pgSz[^>]+w:w=\"11906\"[^>]+w:h=\"16838\"", xml.decode()):
            raise AssertionError("DOCX page size is not explicitly A4.")
    print(f"DOCX validated: {DOCX} ({DOCX.stat().st_size:,} bytes, 9 figures)")


if __name__ == "__main__":
    try:
        validate_docx()
    except Exception as error:
        print(f"DOCX validation failed: {error}", file=sys.stderr)
        raise
