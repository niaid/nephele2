#!/usr/bin/env python3

"""

Format HOMD database for DADA2. See: https://benjjneb.github.io/dada2/training.html

.. code-block:: bash

  Usage: format_homd_dada2.py [-h] [--genus FILE] [--species FILE] fasta taxonomy


positional arguments:
  /fasta
                        FASTA file of sequences from `HOMD download page <http://www.homd.org/index.php?name=seqDownload&file&type=R>`__ .
                        We used *HOMD 16s
                        rRNA RefSeq Version \*\* (Starts from position 9)* file.
                        Aligned FASTA File does not work.
  /taxonomy
                        MOTHUR taxonomy file. Order must match the FASTA file.

optional arguments:
  -h, --help            show this help message and exit
  --genus FILE, -g FILE
                        Output fasta file formatted as DADA2
                        database with genus level header for use with assignTaxonomy. Default will be
                        same as input fasta file with *.train_set.fa* extension.
  --species FILE, -s FILE
                        Output fasta file with species header formatted as DADA2
                        database for use with assignSpecies or addSpecies.
                        Default will be same as input fasta file with
                        *.species_assignment.fa* extension.

"""

import argparse
import re
import sys
import os

genus_ext = ".train_set.fa"
species_ext = ".species_assignment.fa"

def parse_args():
    """Parse arguments"""
    parser = argparse.ArgumentParser(description="Format HOMD database for DADA2. "
                                     "https://benjjneb.github.io/dada2/training.html")
    parser.add_argument("fasta", type=argparse.FileType("r"), metavar = "fasta", help=
                        "FASTA file of sequences from http://www.homd.org/index.php?name=seqDownload&file&type=R . "
                        "We used HOMD 16s rRNA RefSeq Version ** (Starts from position 9). "
                        "Aligned FASTA File does not work.")
    parser.add_argument("taxonomy", type=argparse.FileType("r"), metavar="taxonomy", help=
                        "MOTHUR taxonomy file. Order must match the FASTA file.")
    parser.add_argument("--genus", "-g", type=argparse.FileType("w"), metavar = "FILE", help=
                        "Output fasta file with genus level header formatted as DADA2 database for use with "
                        "assignTaxonomy. Default will be same as input fasta file with " + genus_ext +
                        " extension.")
    parser.add_argument("--species", "-s", type=argparse.FileType("w"), metavar = "FILE", help=
                        "Output fasta file with species level header formatted as DADA2 database for use with "
                        "assignSpecies or addSpecies.  Default will be same as input fasta file with " +
                        species_ext + " extension.")
    return parser.parse_args()


def openoutfiles(fn, ext):
    """Open output fasta files for writing"""
    if fn is None:
        x = re.sub(r"(\.p9)*\.fasta", "", os.path.basename(args.fasta.name)) + ext
        return open(x, "w")
    return fn


if __name__ == "__main__":
    args = parse_args()
    args.genus = openoutfiles(args.genus, genus_ext)
    args.species = openoutfiles(args.species, species_ext)

    with args.fasta, args.taxonomy, args.genus, args.species:
        for fasta_line, tax_line in zip(args.fasta, args.taxonomy):
            ## species file
            sp_header = fasta_line.split(" | ")
            print("{0} {1}".format(sp_header[0], sp_header[1]), file=args.species)
            seq = next(args.fasta)
            print(seq, file=args.species, end="")
            ## genus/whole taxonomy file
            tax_line = tax_line.split()
            if tax_line[0] != sp_header[0][1:]:
                print("{0} {1} fasta ID does not match taxonomy ID.\nExiting."
                      .format(tax_line[0], sp_header[0]), file=sys.stderr)
                sys.exit(1)
            tax = tax_line[1].split(";")
            n6 = min(len(tax) - 1, 6)
            db_header = ">" + ";".join(tax[0:n6]) + ";"
            print(db_header, file=args.genus)
            print(seq, file=args.genus, end="")
