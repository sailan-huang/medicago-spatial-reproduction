from __future__ import annotations

import hashlib
import sys
import urllib.request
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, destination: Path, expected_sha256: str) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and sha256(destination) == expected_sha256.lower():
        print(f"Already present: {destination}")
        return
    temporary = destination.with_suffix(destination.suffix + ".part")
    print(f"Downloading {url}")
    urllib.request.urlretrieve(url, temporary)
    observed = sha256(temporary)
    if observed != expected_sha256.lower():
        temporary.unlink(missing_ok=True)
        raise RuntimeError(
            f"Checksum mismatch for {destination.name}: expected "
            f"{expected_sha256}, observed {observed}"
        )
    temporary.replace(destination)
    print(f"Verified and saved: {destination}")


def main() -> int:
    with (ROOT / "config" / "analysis.yaml").open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)["data"]
    download(
        config["sample_url"],
        ROOT / "data" / "raw" / config["sample_filename"],
        config["sample_sha256"],
    )
    download(
        config["annotation_url"],
        ROOT / "data" / "external" / config["annotation_filename"],
        config["annotation_sha256"],
    )
    download(
        config["supplementary_url"],
        ROOT / "data" / "external" / config["supplementary_filename"],
        config["supplementary_sha256"],
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

