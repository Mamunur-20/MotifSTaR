#!/usr/bin/env python3
"""Download a reproducible subset of public gnomAD REViewer SVGs.

The gnomAD genotype table contains the anonymized ReadvizFilename for each
sample/locus.  This script selects unique primary-locus images, optionally
filters them by PCR protocol and Q score, and downloads the public .svg.gz
objects into one directory per locus.

The selection manifest preserves the available population, genotype, quality,
technical, neurological-case, and public-project metadata for every selected
image so that downstream filtering decisions remain auditable.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import random
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


DEFAULT_BASE_URL = (
    "https://storage.googleapis.com/gnomad-str-public/"
    "release_2024_07/readviz_v2"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download public gnomAD STR REViewer SVG.GZ files"
    )
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--per-locus", type=int, default=100)
    parser.add_argument(
        "--loci",
        default="",
        help="Comma-separated loci; omit to use every locus in the table",
    )
    parser.add_argument(
        "--pcr-protocol",
        choices=("pcr_free", "pcr_plus", "unknown", "any"),
        default="pcr_free",
    )
    parser.add_argument("--min-q", type=float, default=0.0)
    parser.add_argument("--seed", type=int, default=20260901)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--base-url", default=DEFAULT_BASE_URL)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def open_table(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def select_rows(args: argparse.Namespace) -> dict[str, list[dict[str, str]]]:
    requested = {
        value.strip().upper()
        for value in args.loci.split(",")
        if value.strip()
    }
    candidates: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)

    with open_table(args.metadata) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "LocusId",
            "IsAdjacentRepeat",
            "PcrProtocol",
            "Filter",
            "Q",
            "ReadvizFilename",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(
                "Metadata table is missing columns: " + ", ".join(sorted(missing))
            )

        optional_manifest_columns = {
            "IsNeuroCase",
            "ReadLength",
            "PublicProjectId",
            "PublicSampleId",
            "ManualReviewGenotypeQuality",
        }
        missing_optional = optional_manifest_columns.difference(
            reader.fieldnames or []
        )
        if missing_optional:
            print(
                "WARNING: optional metadata columns are absent and will be "
                "blank in the manifest: "
                + ", ".join(sorted(missing_optional))
            )

        for row in reader:
            locus = row["LocusId"].strip()
            filename = row["ReadvizFilename"].strip()
            if not locus or not filename:
                continue
            if requested and locus.upper() not in requested:
                continue
            # Companion components reuse the same SVG filename.  Retaining only
            # the main component prevents duplicate downloads at compound loci.
            if row["IsAdjacentRepeat"].strip().lower() == "true":
                continue
            if row["Filter"].strip().upper() != "PASS":
                continue
            if (
                args.pcr_protocol != "any"
                and row["PcrProtocol"].strip().lower() != args.pcr_protocol
            ):
                continue
            try:
                q_score = float(row["Q"])
            except (TypeError, ValueError):
                continue
            if q_score < args.min_q:
                continue
            candidates[locus][filename] = row

    if requested:
        found = {locus.upper() for locus in candidates}
        missing_loci = sorted(requested.difference(found))
        if missing_loci:
            raise SystemExit(
                "No eligible public images found for: " + ", ".join(missing_loci)
            )

    rng = random.Random(args.seed)
    selected: dict[str, list[dict[str, str]]] = {}
    for locus in sorted(candidates):
        rows = list(candidates[locus].values())
        rng.shuffle(rows)
        selected[locus] = rows[: min(args.per_locus, len(rows))]
    return selected


def image_url(base_url: str, locus: str, filename: str) -> str:
    return "/".join(
        (
            base_url.rstrip("/"),
            urllib.parse.quote(locus, safe=""),
            urllib.parse.quote(filename, safe=""),
        )
    )


def download_one(
    url: str,
    destination: Path,
    retries: int,
) -> tuple[str, str]:
    if destination.exists() and destination.stat().st_size > 0:
        return "existing", str(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".part")
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "gnomAD-STR-public-SVG-downloader/1.0"},
            )
            with urllib.request.urlopen(request, timeout=120) as response:
                with temporary.open("wb") as output:
                    while True:
                        block = response.read(1024 * 1024)
                        if not block:
                            break
                        output.write(block)
            # Google Storage currently serves these .svg.gz objects with
            # transparent decompression.  If a mirror returns real gzip bytes,
            # normalize those too so the local file is always a valid .svg.
            with temporary.open("rb") as downloaded:
                is_gzip = downloaded.read(2) == b"\x1f\x8b"
            if is_gzip:
                normalized = temporary.with_name(temporary.name + ".svg")
                with gzip.open(temporary, "rb") as source, normalized.open("wb") as output:
                    while True:
                        block = source.read(1024 * 1024)
                        if not block:
                            break
                        output.write(block)
                temporary.unlink()
                normalized.replace(destination)
            else:
                temporary.replace(destination)
            return "downloaded", str(destination)
        except (OSError, urllib.error.URLError, urllib.error.HTTPError) as exc:
            last_error = exc
            if temporary.exists():
                temporary.unlink()
            if attempt < retries:
                time.sleep(min(2**attempt, 10))
    return "failed", f"{destination}: {last_error}"


def main() -> int:
    args = parse_args()
    if args.per_locus < 1 or args.workers < 1 or args.retries < 1:
        raise SystemExit("--per-locus, --workers and --retries must be positive")
    if args.min_q < 0:
        raise SystemExit("--min-q must be zero or greater")
    if not args.metadata.is_file():
        raise SystemExit(f"Metadata table does not exist: {args.metadata}")

    selected = select_rows(args)
    args.output.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output / "download_manifest.tsv"
    file_map_path = args.output / "file_map.tsv"
    jobs: list[tuple[str, Path, dict[str, str]]] = []
    with (
        manifest_path.open("w", encoding="utf-8", newline="") as handle,
        file_map_path.open("w", encoding="utf-8", newline="") as map_handle,
    ):
        columns = [
            "sampleID",
            "LocusId",
            "ReadvizFilename",
            "Population",
            "Sex",
            "Age",
            "PcrProtocol",
            "ReadLength",
            "IsNeuroCase",
            "PublicProjectId",
            "PublicSampleId",
            "Genotype",
            "Allele1",
            "Allele2",
            "Filter",
            "Q",
            "ManualReviewGenotypeQuality",
            "url",
            "local_path",
        ]
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        map_writer = csv.DictWriter(
            map_handle,
            fieldnames=("svg_path", "sampleID", "locus", "batch"),
            delimiter="\t",
        )
        map_writer.writeheader()
        for locus, rows in selected.items():
            for row in rows:
                filename = row["ReadvizFilename"].strip()
                sample_id = filename.split(".", 1)[0]
                url = image_url(args.base_url, locus, filename)
                local_filename = filename[:-3] if filename.lower().endswith(".gz") else filename
                destination = args.output / locus / local_filename
                writer.writerow(
                    {
                        **{column: row.get(column, "") for column in columns},
                        "sampleID": sample_id,
                        "LocusId": locus,
                        "ReadvizFilename": filename,
                        "url": url,
                        "local_path": str(destination),
                    }
                )
                map_writer.writerow(
                    {
                        "svg_path": str(destination.relative_to(args.output)),
                        "sampleID": sample_id,
                        "locus": locus,
                        "batch": "gnomAD_population_reference",
                    }
                )
                jobs.append((url, destination, row))

    print(f"Selected {len(jobs):,} unique images across {len(selected)} loci")
    print(f"Selection manifest: {manifest_path}")
    print(f"STR-tool file map: {file_map_path}")
    for locus in sorted(selected):
        print(f"  {locus}: {len(selected[locus])}")
    if args.dry_run:
        print("Dry run: no images downloaded")
        return 0

    counts = defaultdict(int)
    failures: list[str] = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(download_one, url, destination, args.retries): destination
            for url, destination, _ in jobs
        }
        for number, future in enumerate(as_completed(futures), start=1):
            status, message = future.result()
            counts[status] += 1
            if status == "failed":
                failures.append(message)
            if number == 1 or number % 100 == 0 or number == len(futures):
                print(f"Completed {number:,}/{len(futures):,} downloads")

    print(
        "Finished: "
        + ", ".join(f"{key}={value:,}" for key, value in sorted(counts.items()))
    )
    if failures:
        failure_path = args.output / "download_failures.txt"
        failure_path.write_text("\n".join(failures) + "\n", encoding="utf-8")
        print(f"Failures written to: {failure_path}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
