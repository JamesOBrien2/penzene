"""Draws every row of compounds.csv into figures/, one SVG per compound, in the RSC style."""
import csv
import os

import penzene as pz

os.makedirs("figures", exist_ok=True)
with open("compounds.csv", newline="") as f:
    for row in csv.DictReader(f):
        doc = pz.from_smiles(row["smiles"])
        doc.clean()                  # a fresh RDKit layout
        doc.style = "RSC"
        doc.export(f"figures/{row['name']}.svg")
        print(row["name"], doc.formula, round(doc.mw, 2))
