from __future__ import annotations

import hashlib
import sys
import time
from pathlib import Path
from urllib.parse import urlparse

import requests
import yaml

ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(
    url: str,
    destination: Path,
    expected_sha256: str,
    insecure_hosts: set[str],
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and sha256(destination) == expected_sha256.lower():
        print(f"Already present: {destination}")
        return
    temporary = destination.with_suffix(destination.suffix + ".part")
    print(f"Downloading {url}")
    verify_tls = urlparse(url).hostname not in insecure_hosts
    if not verify_tls:
        print(
            "WARNING: TLS verification disabled for a known endpoint; "
            "the downloaded file will still be checked against SHA-256."
        )
    last_error: Exception | None = None
    for attempt in range(1, 7):
        offset = temporary.stat().st_size if temporary.exists() else 0
        headers = {"User-Agent": "medicago-spatial-reproduction/0.1"}
        if offset:
            headers["Range"] = f"bytes={offset}-"
            print(f"Resuming at {offset:,} bytes (attempt {attempt}/6)")
        try:
            with requests.get(
                url,
                stream=True,
                timeout=(30, 180),
                headers=headers,
                verify=verify_tls,
            ) as response:
                if offset and response.status_code == 416:
                    break
                response.raise_for_status()
                mode = "ab" if offset and response.status_code == 206 else "wb"
                with temporary.open(mode) as handle:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            handle.write(chunk)
        except (requests.RequestException, OSError) as error:
            last_error = error
            print(f"Download interrupted: {error}")
            time.sleep(min(30, attempt * 3))
            continue
        if temporary.exists() and sha256(temporary) == expected_sha256.lower():
            break
        last_error = RuntimeError("Downloaded file does not match the pinned digest")
        time.sleep(min(30, attempt * 3))
    else:
        raise RuntimeError(f"Could not download {url}") from last_error

    observed = sha256(temporary) if temporary.exists() else ""
    if observed != expected_sha256.lower():
        raise RuntimeError(
            f"Checksum mismatch for {destination.name}: expected "
            f"{expected_sha256}, observed {observed}"
        )
    temporary.replace(destination)
    print(f"Verified and saved: {destination}")


def main() -> int:
    with (ROOT / "config" / "analysis.yaml").open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)["data"]
    insecure_hosts = set(config.get("allow_insecure_tls_hosts", []))
    download(
        config["sample_url"],
        ROOT / "data" / "raw" / config["sample_filename"],
        config["sample_sha256"],
        insecure_hosts,
    )
    download(
        config["annotation_url"],
        ROOT / "data" / "external" / config["annotation_filename"],
        config["annotation_sha256"],
        insecure_hosts,
    )
    download(
        config["supplementary_url"],
        ROOT / "data" / "external" / config["supplementary_filename"],
        config["supplementary_sha256"],
        insecure_hosts,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
