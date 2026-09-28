# NMR shift table

`hose.tar.xz` holds `hose.tsv`: 13C and 1H chemical shifts keyed by HOSE code, which Penzene
uses to predict shifts (View → Predicted NMR Shifts, #403).

Contains information from nmrshiftdb2 (www.nmrshiftdb.org), which is made available here under
the nmrshiftdb2 Database License (https://nmrshiftdb.nmr.uni-koeln.de/nmrshiftdbhtml/nmrshiftdb2datalicense.txt,
a copy in [LICENSE](LICENSE)). The table is a Derivative Database of nmrshiftdb2 and is licensed
under that licence, not Penzene's GPL-3.

How it's made (§4.7): `cmake/nmr-table.py` reads nmrshiftdb2's public export
(nmrshiftdb2withsignals.sd, from https://sourceforge.net/projects/nmrshiftdb2/files/data/) and
records, for every assigned 13C and 1H shift, the HOSE code of its atom at 1 to 4 spheres, as
Penzene computes them. Each line is `<spheres> <code>`, then the median 13C shift and its count,
then the median 1H shift (the atom's H, averaged per spectrum) and its count. Run the script again
to rebuild it from a newer export.

This table: nmrshiftdb2withsignals.sd downloaded 2026-09-28, SHA-256
0e86688360e23c88ccf0eb82a1251315fa57ec6f5b376dc8886f3311793a0afe, 43,502 compounds with 13C or
1H spectra, 363,296 codes (the same notice and provenance head hose.tsv).
