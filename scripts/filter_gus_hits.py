"""Strict project thresholds for the 15-column InvestiGUT DIAMOND format."""
import argparse
import csv
import math
from collections import Counter


def qualifies(fields):
    if len(fields) != 15:
        raise ValueError(f"Expected 15 fields, found {len(fields)}")
    identity, qcov, scov, evalue = map(float, (fields[2], fields[12], fields[13], fields[10]))
    if not all(math.isfinite(x) for x in (identity, qcov, scov, evalue)):
        raise ValueError("Non-finite alignment metric")
    if not all(0 <= x <= 100 for x in (identity, qcov, scov)) or evalue < 0:
        raise ValueError("Alignment metric outside valid range")
    return identity >= 95 and qcov >= 90 and scov >= 90 and evalue < 1e-5


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("query", "input", "output", "qc"):
        parser.add_argument(f"--{name}", required=True)
    args = parser.parse_args()
    with open(args.query) as handle:
        ids = [line[1:].split()[0] for line in handle if line.startswith(">")]
    if len(ids) != 11 or len(set(ids)) != 11:
        raise ValueError("Expected 11 unique GUS query IDs")
    raw, kept = Counter(), Counter()
    with open(args.input) as source, open(args.output, "w") as target:
        for n, line in enumerate(source, 1):
            fields = line.rstrip("\n").split("\t")
            if fields[0] not in ids:
                raise ValueError(f"Line {n}: unknown query {fields[0]}")
            raw[fields[0]] += 1
            if qualifies(fields):
                target.write(line)
                kept[fields[0]] += 1
    with open(args.qc, "w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["query_id", "raw_hit_rows", "qualified_hit_rows", "excluded_hit_rows"])
        for query in ids:
            writer.writerow([query, raw[query], kept[query], raw[query] - kept[query]])


if __name__ == "__main__":
    main()
