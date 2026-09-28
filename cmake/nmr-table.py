"""Build resources/nmr/hose.tar.xz, the NMR shift table behind View → Predicted NMR Shifts (#403).

From nmrshiftdb2's public export (nmrshiftdb2withsignals.sd, nmrshiftdb2 Database License): every
assigned 13C and 1H shift, keyed by the HOSE code of its atom (the carbon, or the atom the H is on)
at 1 to 4 spheres. Each line of hose.tsv is "<spheres> <code>\t13C\tcount\t1H\tcount": the median
shift in ppm and how many shifts it's from (empty without any), sorted by key for binary search.

The codes come from Penzene's own HOSE code function (chem::hoseCodes, through the Python module),
so the app's lookups match the table exactly. Run from the repository root after a build:

    PYTHONPATH=build/python python cmake/nmr-table.py [nmrshiftdb2withsignals.sd] [--evaluate]

Without a file it downloads the current export into build/. --evaluate also holds out one
compound in ten (by InChIKey) and prints the mean absolute error of predicting them.
"""
import hashlib
import io
import re
import statistics
import sys
import tarfile
import urllib.request
from collections import defaultdict
from pathlib import Path

from penzene import _penzene

URL = "https://sourceforge.net/projects/nmrshiftdb2/files/data/nmrshiftdb2withsignals.sd/download"
ROOT = Path(__file__).resolve().parent.parent
SPHERES = 4
FIELD = re.compile(r"^> +<([^>]+)>\n(.*?)(?:\n\n|\n?\Z)", re.M | re.S)


def records(path):
    """(InChIKey, {atom: [13C]}, {atom: [1H]}, codes) per compound with a 13C or 1H spectrum."""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    skipped = 0
    for rec in text.split("$$$$\n"):
        block, _, tail = rec.partition("M  END")
        if not tail:
            continue
        fields = {k: v for k, v in FIELD.findall(tail)}
        shifts = {"13C": defaultdict(list), "1H": defaultdict(list)}
        for name, value in fields.items():
            m = re.match(r"Spectrum (13C|1H) \d+$", name)
            if not m:
                continue
            per_atom = defaultdict(list)  # one value per atom per spectrum: a CH3's three H count once
            for signal in value.strip().split("|"):
                parts = signal.split(";")
                if len(parts) == 3 and parts[2].strip().isdigit():
                    per_atom[int(parts[2])].append(float(parts[0]))
            for atom, values in per_atom.items():
                shifts[m.group(1)][atom].append(sum(values) / len(values))
        if not shifts["13C"] and not shifts["1H"]:
            continue
        codes = _penzene._hose_codes(block + "M  END\n", SPHERES)
        if not codes:
            skipped += 1
            continue
        for atom in [a for a in shifts["1H"] if a < len(codes) and not codes[a]]:  # on a drawn H: its atom's
            if (on := hydrogen_partner(block, atom)) is not None:
                shifts["1H"][on] += shifts["1H"].pop(atom)
        yield fields.get("INChI key", "").strip(), shifts["13C"], shifts["1H"], codes
    print(f"skipped {skipped} records RDKit can't read", file=sys.stderr)


def hydrogen_partner(block, atom):
    """The atom an H atom (0-based) is bonded to, from a V2000 block; None if it isn't bonded once."""
    lines = block.split("\n")
    counts = next(i for i, line in enumerate(lines) if "V2000" in line)
    atoms, bonds = int(lines[counts][:3]), int(lines[counts][3:6])
    partners = [int(b[3:6]) - 1 if int(b[:3]) - 1 == atom else int(b[:3]) - 1
                for b in lines[counts + 1 + atoms:counts + 1 + atoms + bonds] if atom + 1 in (int(b[:3]), int(b[3:6]))]
    return partners[0] if len(partners) == 1 else None


def build(data):
    """{"<spheres> <code>": ([13C], [1H])} from records."""
    table = defaultdict(lambda: ([], []))
    bad = defaultdict(int)
    for _, c13, h1, codes in data:
        for column, shifts in enumerate((c13, h1)):
            for atom, values in shifts.items():
                if atom >= len(codes) or not codes[atom] or (column == 0 and not codes[atom][0].startswith("C-")):
                    bad[column] += 1
                    continue
                for n, code in enumerate(codes[atom], 1):
                    table[f"{n} {code}"][column].extend(values)
    print(f"skipped assignments to no atom, an H atom, or (13C) not a carbon: 13C {bad[0]}, 1H {bad[1]}",
          file=sys.stderr)
    return table


def predict(table, codes, column):
    """(median shift, spheres matched), from the most spheres the table has."""
    for n in range(len(codes), 0, -1):
        values = table.get(f"{n} {codes[n - 1]}", ([], []))[column]
        if values:
            return statistics.median(values), n
    return None, 0


def evaluate(data):
    held = lambda key: int(hashlib.md5(key.split("-")[0].encode()).hexdigest(), 16) % 10 == 0
    table = build(r for r in data if not held(r[0]))
    for column, name in enumerate(("13C", "1H")):
        errors = defaultdict(list)
        missed = 0
        for key, c13, h1, codes in data:
            if not held(key):
                continue
            for atom, values in (c13, h1)[column].items():
                if atom >= len(codes) or not codes[atom]:
                    continue
                shift, n = predict(table, codes[atom], column)
                if n:
                    errors[n] += [abs(shift - v) for v in values]
                else:
                    missed += len(values)
        every = [e for n in errors for e in errors[n]]
        print(f"{name}: MAE {statistics.mean(every):.2f} ppm over {len(every)} held-out shifts, {missed} unpredicted")
        for n in sorted(errors, reverse=True):
            print(f"  {n} spheres: MAE {statistics.mean(errors[n]):.2f} ppm, {len(errors[n])} shifts")


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    sd = Path(args[0]) if args else ROOT / "build" / "nmrshiftdb2withsignals.sd"
    if not sd.exists():
        print("downloading", URL, file=sys.stderr)
        urllib.request.urlretrieve(URL, sd)
    data = list(records(sd))
    if "--evaluate" in sys.argv:
        evaluate(data)
    fmt = lambda values: (f"{statistics.median(values):.2f}", str(len(values))) if values else ("", "")
    sha = hashlib.sha256(sd.read_bytes()).hexdigest()
    lines = ["# Contains information from nmrshiftdb2 (www.nmrshiftdb.org), which is made available here under the "
             "nmrshiftdb2 Database License (https://nmrshiftdb.nmr.uni-koeln.de/nmrshiftdbhtml/nmrshiftdb2datalicense.txt).",
             f"# Built by Penzene's cmake/nmr-table.py from {sd.name} (SHA-256 {sha}, {len(data)} compounds with 13C or 1H)."]
    table = build(data)
    lines += ["\t".join((key, *fmt(table[key][0]), *fmt(table[key][1]))) for key in sorted(table)]
    tsv = ("\n".join(lines) + "\n").encode()
    out = ROOT / "resources" / "nmr" / "hose.tar.xz"
    out.parent.mkdir(exist_ok=True)
    with tarfile.open(out, "w:xz", format=tarfile.USTAR_FORMAT) as tar:  # reproducible: no times or owners
        info = tarfile.TarInfo("hose.tsv")
        info.size, info.mode = len(tsv), 0o644
        tar.addfile(info, io.BytesIO(tsv))
    print(f"wrote {out.relative_to(ROOT)}: {len(table)} codes, {out.stat().st_size / 1e6:.1f} MB", file=sys.stderr)


if __name__ == "__main__":
    main()
