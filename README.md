# Protein Research
Summer Research at the National Institute of Standards and Technology (NIST), Summer 2016

## Background
Glycosylation is the process in which a carbohydrate attaches to a protein or lipid, and it is essential to the folding and function of proteins. Therefore, defects in glycosylation can cause severe diseases and even cancers. Proteins consist of long strings of amino acids, and our research focused mainly on N-linked glycosylation sites, which have an N-X-S/T/C motif (where X is not Proline). For N-linked glycosylation, not all potential glycosylation sites become glycosylated, so this program was designed to predict which glycosylation sites are actually important based on how well the N-X-S/T/C motif is conserved across several species.

## Project Details
This program data-mines the UniProt website for sequences that have experimentally proven n-linked glycosylation sites. It then expands the search and aligns all sequences of a given protein using the Clustal Omega alignment software. From there, it performs statistical analyses on the alignments to predict significant potential glycosylation sites. The data can also be used to generate other statistics, like amino acid (X) frequency, alignment of flanking regions, and mutation rates (i.e. from N to D).

The program performs several optimizations in order to quickly compile and sort large quantities of data, such as specifying which glycosylation sites already have experimentally proven data, filtering out duplicate organisms to find the "most representative" sequence, and distinguishing between N-X-S, N-X-T, and N-X-C sites. It also stores all the information retreived in an easily retrievable format to allow for quick and easy analysis.

1). info_tester.py queries the UniProt database to data-mine and parse all the sequencing data. It also generates gene tables and fasta files for the sequence information.

2). alignment_tester.py runs the Clustal Omega alignment softare, and produces a clean, "pretty printed" alignment for all the sequences of the proteins.

3). analysis_tester.py performs the statistical analyses on the sequencing data. While it is not necessarily an executable file, it has separate function calls from analysis.py (where the statistical functions are actually written) and a rough pipeline structure.

## Data Sources
Protein sequencing data was taken from [UniProt](http://www.uniprot.org/).

## Technology
This program was written using Python 3. It uses the BioPython, tqdm, and Bioservices packages, and the Clustal Omega alignment software. In order to generate colored html files, the [ansifilter](http://www.andre-simon.de/doku/ansifilter/en/ansifilter.php) program is run from the command line in order to generate the visual html alignment files.
