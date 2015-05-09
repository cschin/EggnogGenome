EggnogGenome uses best BLAST hit and the Eggnog database to annotate functions on a query set of proteins.
Eggnog uses the `COG functions notation <http://eggnogdb.embl.de/download/latest/eggnog4.functional_categories.txt>`_

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