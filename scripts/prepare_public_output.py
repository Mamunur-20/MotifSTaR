#!/usr/bin/env python3

import argparse
import gzip
import shutil
import zipfile
from pathlib import Path


REPLACEMENTS = (
    (
        "/g/data/kp73/mr3328/MND_project/STR_tool/data/"
        "gene_strands_companion_locus.xlsx",
        "config/gene_strands_companion_locus.xlsx",
    ),
    (
        "/g/data/kp73/mr3328/MND_project/STR_tool/data/"
        "variant_catalog_without_offtargets.GRCh38.json",
        "config/variant_catalog_without_offtargets.GRCh38.json",
    ),
    (
        "/g/data/kp73/mr3328/MND_project/STR_tool/data/"
        "gnomAD_100_per_locus/file_map.tsv",
        "external_inputs/gnomAD_100_per_locus/file_map.tsv",
    ),
    (
        "/g/data/kp73/mr3328/MND_project/STR_tool/data/"
        "gnomAD_100_per_locus/",
        "external_inputs/gnomAD_100_per_locus/",
    ),
    (
        "/g/data/kp73/mr3328/MND_project/STR_tool/data/"
        "gnomAD_100_per_locus",
        "external_inputs/gnomAD_100_per_locus",
    ),
)

PLAIN_SUFFIXES = {".csv", ".tsv", ".json", ".log", ".txt", ".md"}
PRIVATE_MARKERS = ("/g/data/", "mr3328")


def replace_paths(text):
    replacements = 0
    for old, new in REPLACEMENTS:
        replacements += text.count(old)
        text = text.replace(old, new)
    return text, replacements


def clean_plain_file(path):
    text = path.read_text(encoding="utf-8")
    updated, count = replace_paths(text)
    if count:
        path.write_text(updated, encoding="utf-8", newline="")
    return count


def clean_gzip_file(path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        text = handle.read()

    updated, count = replace_paths(text)
    if not count:
        return 0

    temporary = path.with_name(path.name + ".tmp")
    with gzip.open(temporary, "wt", encoding="utf-8", newline="") as handle:
        handle.write(updated)

    temporary.replace(path)
    return count


def clean_xlsx_file(path):
    temporary = path.with_name(path.name + ".tmp")
    replacements = 0

    with zipfile.ZipFile(path, "r") as source:
        with zipfile.ZipFile(
            temporary, "w", compression=zipfile.ZIP_DEFLATED
        ) as destination:
            for member in source.infolist():
                data = source.read(member.filename)

                if member.filename.endswith((".xml", ".rels")):
                    try:
                        text = data.decode("utf-8")
                    except UnicodeDecodeError:
                        pass
                    else:
                        text, count = replace_paths(text)
                        replacements += count
                        data = text.encode("utf-8")

                destination.writestr(member, data)

    temporary.replace(path)
    return replacements


def contains_private_marker(path):
    if path.suffix.lower() in PLAIN_SUFFIXES:
        text = path.read_text(encoding="utf-8", errors="ignore")
        return any(marker in text for marker in PRIVATE_MARKERS)

    if path.name.lower().endswith(".tsv.gz"):
        with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as handle:
            text = handle.read()
        return any(marker in text for marker in PRIVATE_MARKERS)

    if path.suffix.lower() == ".xlsx":
        with zipfile.ZipFile(path, "r") as workbook:
            for name in workbook.namelist():
                if not name.endswith((".xml", ".rels")):
                    continue
                text = workbook.read(name).decode("utf-8", errors="ignore")
                if any(marker in text for marker in PRIVATE_MARKERS):
                    return True

    return False


def main():
    parser = argparse.ArgumentParser(
        description="Prepare a portable public copy of MotifSTaR output."
    )
    parser.add_argument("--source", required=True)
    parser.add_argument("--destination", required=True)
    args = parser.parse_args()

    source = Path(args.source).resolve()
    destination = Path(args.destination).resolve()

    if not source.is_dir():
        raise SystemExit(f"Source directory does not exist: {source}")

    if destination.exists():
        raise SystemExit(
            f"Destination already exists: {destination}\n"
            "Use a new destination name or remove the incomplete copy."
        )

    shutil.copytree(source, destination)

    changed_files = 0
    replacement_count = 0

    for path in destination.rglob("*"):
        if not path.is_file():
            continue

        count = 0

        if path.name.lower().endswith(".tsv.gz"):
            count = clean_gzip_file(path)
        elif path.suffix.lower() == ".xlsx":
            count = clean_xlsx_file(path)
        elif path.suffix.lower() in PLAIN_SUFFIXES:
            count = clean_plain_file(path)

        if count:
            changed_files += 1
            replacement_count += count
            print(f"Cleaned {count} path occurrence(s): {path}")

    unsafe_files = [
        path
        for path in destination.rglob("*")
        if path.is_file() and contains_private_marker(path)
    ]

    if unsafe_files:
        print("\nPrivate paths remain in:")
        for path in unsafe_files:
            print(path)
        raise SystemExit("Public-output preparation failed.")

    print()
    print(f"Public output created: {destination}")
    print(f"Files changed: {changed_files}")
    print(f"Path occurrences replaced: {replacement_count}")
    print("No private Gadi paths remain.")


if __name__ == "__main__":
    main()