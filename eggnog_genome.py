#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EggnogGenome uses best BLAST hit and the Eggnog database to annotate functions on a query set of proteins.
Eggnog uses the `COG functions notation <http://eggnogdb.embl.de/download/latest/eggnog4.functional_categories.txt>`

Expected workflow:
1. Download Eggnog database. Required files are:
    * eggnog4.proteins.core_periphery.fa.gz: protein FASTA that will be used as BLAST database
    * data/NOG/NOG.members.tsv.gz: table mapping Eggnog proteins to their COG function(s)

2. Transform the Eggnog proteins FASTA into BLAST database.
   Eg: makeblastdb -dbtype prot -in eggnog4.proteins.core_periphery.fa -out eggnog.blastdb

3. BLASTp your query proteins against Eggnog BLAST database. Eg:
   * Legacy BLAST: blastall -p blastp -d eggnog.blastdb -i query_proteins.fa -o matches.blast
     -e .01 -m 7 -a 60 -K 1
   * BLAST+: blastp -db eggnog.blastdb -query query_proteins.fa -out matches.blast -evalue .01
     -max_hsps 1 -outfmt 5 -num_alignments 1 -num_threads 60

4. Run eggnog_genome.py


`Eggnog homepage <http://eggnogdb.embl.de/#/app/home>`_
`Eggnog downloads <http://eggnogdb.embl.de/download/latest/>`_
`NCBI BLAST <http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`_
"""

__license = """GPL v3

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
__version__ = "0.1"
__author__ = "Sébastien Gélis-Jeanvoine"
__copyright__ = "Copyright 2015, Sébastien Gélis-Jeanvoine"
__credits__ = ["Sébastien Gélis-Jeanvoine"]
__maintainer__ = "Sébastien Gélis-Jeanvoine"
__email__ = "sebastien@gelis.ch"

import argparse
from collections import Counter, defaultdict
import re
import sys
from Bio.Blast import NCBIXML


def parse_eggnog(members_file):
    """
    Read an Eggnog members file (typically NOG.members.tsv) and return a
    dictionary mapping Eggnog proteins to functions.

    :param members_file: path to file containing Eggnog members
    :type members_file: str

    :return: dictionary mapping Eggnog proteins to functions
    :rtype: defaultdict
    """
    print("Parsing Eggnog members...")

    prots_to_funcs = defaultdict(str)

    with open(members_file, "r") as members_handle:
        for line in members_handle:
            splitline = line.split("\t")
            func, prots = splitline[4], splitline[5].rstrip().split(",")
            for prot in prots:
                prots_to_funcs[prot] = func

    return prots_to_funcs


def parse_blast(blast_file, prots_to_funcs):
    """
    Reads a BLAST XML output and returns a dictionary mapping hits to
    functions.

    :param blast_file: path to BLAST XML file
    :type blast_file: str
    :param prots_to_funcs: dictionary mapping Eggnog proteins to functions
    :type prots_to_funcs: defaultdict

    :return: dictionary mapping BLAST hits to Eggnog functions
    :rtype: dict
    """
    print("Parsing BLAST matches...")

    gene_name = re.compile("\[gene=(.*?)\]")
    gene = lambda x: gene_name.findall(x)[0]

    hits_to_funcs = {}

    try:
        blast_filehandle = open(blast_file, "r")
        blast_handle = NCBIXML.parse(blast_filehandle)
    except FileNotFoundError:
        sys.stderr.write("Error: could not find file {}\n".format(blast_file))
    except ValueError:
        sys.stderr.write("Error: invalid XML BLAST file: {}".format(blast_file))
    else:

        for record in blast_handle:
            try:
                aln = record.alignments[0]
            except IndexError:
                hits_to_funcs[gene(record.query)] = ""
            else:
                sbjct = aln.title.split()[1]
                hits_to_funcs[gene(record.query)] = prots_to_funcs[sbjct]

        blast_filehandle.close()

    return hits_to_funcs


def log(hits_to_funcs):
    """
    Write some metrics to a log file for user's information.

    :param hits_to_funcs: dictionary mapping BLAST hits to Eggnog functions
    :type hits_to_funcs: dict

    :return:
    :rtype: None
    """
    print("Writing log file...")

    funcs_counter = Counter()
    hypothetical = funcs_counter["S"] + funcs_counter["R"]
    nonempty = len(hits_to_funcs) - funcs_counter[""]

    with open("eggnog.log", "w") as outf:
        outf.write("Total query proteins: {}\n".format(len(hits_to_funcs)))
        outf.write("Annotated: {0} ({1:0.2f}%)\n".format(nonempty, nonempty / len(hits_to_funcs) * 100))
        outf.write("Hypothetical: {0} ({1:0.2f}%)\n".format(hypothetical, hypothetical / len(hits_to_funcs) * 100))


def histogram(hits_to_funcs):
    """
    Write the functions histogram. Can be used subsequently to pie chart genome.

    :param hits_to_funcs: dictionary mapping BLAST hits to Eggnog functions
    :type hits_to_funcs: dict

    :return:
    :rtype: None
    """
    print("Writing histogram...")

    funcs_histogram = Counter(hits_to_funcs.values())
    with open("eggnog_histogram.csv", "w") as outf:
        for func in funcs_histogram.elements():
            outf.write("{0};{1}\n".format(func, funcs_histogram[func]))


def main():
    """
    Parse args and start analyses with helper functions.
    """
    # Parse command line arguments
    argparser = argparse.ArgumentParser(prog="eggnog_genome.py",
                                        description="Eggnog your genome")
    argparser.add_argument("-v", "--version", action="version",
                           version="%(prog)s 0.1")
    argparser.add_argument("-b", "--blast", nargs=1, action="store",
                           dest="blast", required=True, help="Input BLAST results (in XML format)")
    argparser.add_argument("-m", "--members", nargs=1, action="store",
                           dest="members_file", required=True, help="Eggnog members table")
    argparser.add_argument("-o", "--output", nargs=1, action="store", dest="output",
                           default="eggnog_genome_output.csv", help="Filename to write result table to. "
                                                                    "Default: eggnog_genome_output.csv")
    argparser.add_argument("-s", "--separator", nargs=1, action="store", dest="sep",
                           choices=[";", "tab"], default=";", help="Columns separator in output file. "
                                                                   "Default: ;")
    argparser.add_argument("-l", "--log", action="store_true", dest="log_flag", help="Log some metrics")
    argparser.add_argument("--hist", action="store_true", dest="hist_flag", help="Create functions histogram")
    args = argparser.parse_args()

    # Start analyses
    prots_to_funcs = parse_eggnog(args.members_file[0])
    hits_to_funcs = parse_blast(args.blast[0], prots_to_funcs)
    if args.log_flag:
        log(hits_to_funcs)
    if args.hist_flag:
        histogram(hits_to_funcs)

    # Write final table to output file
    print("Writing results to output file...")
    with open(args.output, "w") as outf:
        for hit in hits_to_funcs.keys():
            outf.write("{0}{1}{2}\n".format(hit,
                                            args.sep[0] if args.sep[0] != "tab" else "\t",
                                            hits_to_funcs[hit]))
    print("All done!")

if __name__ == "__main__":
    main()
