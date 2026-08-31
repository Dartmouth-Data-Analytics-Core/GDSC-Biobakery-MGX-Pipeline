#!/usr/bin/env python3
"""
stratify_abundance.py

Pull a single taxonomic rank out of a merged MetaPhlAn abundance table and write
it back out with real sample names as the column header.

This replaces the grep/sed chain that used to live in rule extract_genus_species.
That chain dropped the header line, because the header contains neither "s__"
nor "g__", so the stratified tables came out with no sample names at all and the
R scripts downstream had to guess the column order from the sample sheet. If the
merge order and the sample sheet order ever disagreed, every sample was silently
relabelled and nothing errored.

Here the sample sheet is used to check the names rather than to supply them. If
the columns in the table and the sample_ids in the sheet do not match exactly,
this exits non-zero instead of guessing.

Standard library only.

Input format, as written by merge_metaphlan_tables.py (MetaPhlAn 4.0.6):

    #mpa_vOct22_CHOCOPhlAnSGB_202212
    clade_name	S1-profiled_metagenome	S2-profiled_metagenome
    k__Bacteria	100.0	100.0
    k__Bacteria|p__Firmicutes	60.1	55.4
    ...

Note that the merge keeps only the relative_abundance column from each profile
(it reads usecols=[0,2]), so the extra columns produced by
metaphlan_analysis_type: "rel_ab_w_read_stats" do not reach this point. Values
here are percentages, 0 to 100, not fractions.
"""

import argparse
import csv
import os
import sys

#----- MetaPhlAn rank prefixes, shallowest to deepest
RANK_PREFIXES = {
    "kingdom": "k__",
    "phylum": "p__",
    "class": "c__",
    "order": "o__",
    "family": "f__",
    "genus": "g__",
    "species": "s__",
    "sgb": "t__",
}

#----- Suffix metaphlan/merge leaves on each column name.
# rule taxonomy writes {sample}-profiled_metagenome.txt, and
# merge_metaphlan_tables.py names each column after the file basename.
DEFAULT_SUFFIX = "-profiled_metagenome"


def die(message):
    sys.stderr.write("stratify_abundance.py: error: {}\n".format(message))
    sys.exit(1)


def log(message):
    sys.stderr.write("stratify_abundance.py: {}\n".format(message))


def read_merged_table(path):
    """Return (comments, header_fields, data_rows) from a merged MetaPhlAn table.

    Leading '#' lines are treated as comments. The first line that is not a
    comment is the column header.
    """
    comments = []
    header = None
    rows = []

    with open(path) as handle:
        for lineno, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")

            if not line.strip():
                continue

            if header is None and line.startswith("#"):
                comments.append(line)
                continue

            fields = line.split("\t")

            if header is None:
                header = fields
                continue

            if len(fields) != len(header):
                die(
                    "{}: line {} has {} fields, header has {}".format(
                        path, lineno, len(fields), len(header)
                    )
                )

            rows.append(fields)

    if header is None:
        die("{}: no header line found, file is empty or all comments".format(path))

    if len(header) < 2:
        die(
            "{}: header has {} column(s), expected clade_name plus at least one "
            "sample".format(path, len(header))
        )

    return comments, header, rows


def clean_sample_names(columns, suffix):
    """Strip the profile-file suffix off each column name."""
    cleaned = []
    for column in columns:
        name = column.strip()
        if suffix and name.endswith(suffix):
            name = name[: -len(suffix)]
        cleaned.append(name)
    return cleaned


def read_sample_ids(path):
    """Read the sample_id column from the pipeline sample sheet, in file order."""
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)

        if reader.fieldnames is None:
            die("{}: sample sheet is empty".format(path))

        fieldnames = [f.strip() for f in reader.fieldnames]
        if "sample_id" not in fieldnames:
            die(
                "{}: no 'sample_id' column, found: {}".format(
                    path, ", ".join(fieldnames)
                )
            )

        sample_ids = []
        for row in reader:
            value = (row.get("sample_id") or "").strip()
            if value:
                sample_ids.append(value)

    if not sample_ids:
        die("{}: no sample rows found".format(path))

    duplicates = sorted({s for s in sample_ids if sample_ids.count(s) > 1})
    if duplicates:
        die("{}: duplicate sample_id values: {}".format(path, ", ".join(duplicates)))

    return sample_ids


def resolve_column_order(table_samples, sample_ids, table_path, sheet_path):
    """Check the table columns against the sample sheet, return sheet-order indices.

    Both directions are checked, so a sample missing from the table and a stray
    extra column are both caught.
    """
    duplicates = sorted({s for s in table_samples if table_samples.count(s) > 1})
    if duplicates:
        die(
            "{}: duplicate sample columns after cleaning names: {}".format(
                table_path, ", ".join(duplicates)
            )
        )

    in_table = set(table_samples)
    in_sheet = set(sample_ids)

    missing = [s for s in sample_ids if s not in in_table]
    extra = [s for s in table_samples if s not in in_sheet]

    if missing or extra:
        message = "sample names in {} do not match {}".format(table_path, sheet_path)
        if missing:
            message += "\n  in sample sheet but not in table: {}".format(
                ", ".join(missing)
            )
        if extra:
            message += "\n  in table but not in sample sheet: {}".format(
                ", ".join(extra)
            )
        message += "\n  table columns:  {}".format(", ".join(table_samples))
        message += "\n  sample sheet:   {}".format(", ".join(sample_ids))
        die(message)

    position = {name: i for i, name in enumerate(table_samples)}
    return [position[s] for s in sample_ids]


def rows_at_rank(rows, prefix):
    """Keep rows whose deepest lineage element sits at the requested rank.

    MetaPhlAn writes one row per level with the full pipe-joined lineage, so the
    prefix on the last element identifies that row's rank exactly. Testing the
    last element is what keeps deeper levels out: asking for species will not
    pick up t__ SGB rows, and asking for genus will not pick up s__ rows.
    """
    kept = []
    for fields in rows:
        clade = fields[0].strip()
        leaf = clade.split("|")[-1]
        if leaf.startswith(prefix):
            kept.append((leaf, fields[1:]))
    return kept


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract one taxonomic rank from a merged MetaPhlAn abundance table, "
            "keeping sample names on the columns."
        )
    )
    parser.add_argument(
        "--input", required=True, help="merged MetaPhlAn abundance table"
    )
    parser.add_argument("--output", required=True, help="stratified table to write")
    parser.add_argument(
        "--rank",
        required=True,
        choices=sorted(RANK_PREFIXES),
        help="taxonomic rank to extract",
    )
    parser.add_argument(
        "--samples",
        help=(
            "pipeline sample sheet. When given, the table columns are checked "
            "against its sample_id column and written in sample sheet order."
        ),
    )
    parser.add_argument(
        "--sample-suffix",
        default=DEFAULT_SUFFIX,
        help="suffix to strip off column names (default: %(default)s)",
    )
    parser.add_argument(
        "--index-name",
        default="clade_name",
        help="name for the first output column (default: %(default)s)",
    )
    args = parser.parse_args()

    prefix = RANK_PREFIXES[args.rank]

    comments, header, rows = read_merged_table(args.input)

    for comment in comments:
        if comment.startswith("#mpa_"):
            log("profile database: {}".format(comment.lstrip("#")))

    table_samples = clean_sample_names(header[1:], args.sample_suffix)

    if args.samples:
        sample_ids = read_sample_ids(args.samples)
        order = resolve_column_order(
            table_samples, sample_ids, args.input, args.samples
        )
        out_samples = sample_ids
    else:
        log("no --samples given, column names taken from the table as-is")
        order = list(range(len(table_samples)))
        out_samples = table_samples

    kept = rows_at_rank(rows, prefix)

    if not kept:
        die(
            "no rows at rank '{}' ({}) in {}. Check that the profiles were built "
            "against the expected database.".format(args.rank, prefix, args.input)
        )

    seen = {}
    for leaf, _ in kept:
        seen[leaf] = seen.get(leaf, 0) + 1
    repeated = sorted(name for name, count in seen.items() if count > 1)
    if repeated:
        log(
            "warning: {} clade name(s) appear more than once at this rank, "
            "for example {}".format(len(repeated), repeated[0])
        )

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.output, "w") as handle:
        handle.write("\t".join([args.index_name] + out_samples) + "\n")
        for leaf, values in kept:
            ordered = [values[i] for i in order]
            handle.write("\t".join([leaf] + ordered) + "\n")

    log(
        "{}: {} rows at rank '{}' across {} samples -> {}".format(
            os.path.basename(args.input),
            len(kept),
            args.rank,
            len(out_samples),
            args.output,
        )
    )


if __name__ == "__main__":
    main()
