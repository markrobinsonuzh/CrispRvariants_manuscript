# CrispRvariants_manuscript

This repository contains the code used in 
"CrispRVariants: precisely charting the mutation spectrum in genome engineering experiments"
by [Lindsay *et al*](http://biorxiv.org/content/early/2015/12/10/034140).

This repository also contains the Sanger sequencing data from Burger *et al* referred
to in the paper.  Other data used in the paper is publically available.  See the pdf
for further details of where the data can be assessed and how these scripts are used.

## Additional scripts

In addition to several freely available software packages described in the pdf document,
this code relies on the following scripts being in the src directory:

+[psl2sam.pl] (https://github.com/samtools/samtools/blob/develop/misc/psl2sam.pl) is part
of the [samtools](https://github.com/samtools/) suite.
+samGetSEQfast.pl from the [rackJ](http://rackj.sourceforge.net/) toolkit. 

We also include a script 
+ ampliconDIV_minimal.sh in the src directory, which we run from the ampliconDIVder 
directory.  This code is an excerpt from the ampliconDIVider_driver.sh, 
written by Matthew LaFave.

## Annotation

### Files for AmpliconDIVider analysis

### Files for CRISPResso analysis

+ CRISPResso_pooled_layout.txt contains the commands used for running CRISPRessoPooled.  The
matching script is src/run_crispresso_pooled.sh
+ shah_custom_amplicons.txt contains the custom amplicon sequences inferred from the differences
from the mapped sequences to the danRer7 reference, which was used for running CRISPResso.

### Metadata files

The following files refer to data from "Rapid reverse genetic screening using CRISPR in zebrafish"
by Shah *et al*:

+shah_guides.bed - the locations of the guides in the danRer7/Zv9 genome.
+shah_primers.bed - the range spanned by the PCR primers, i.e. the amplicon ranges
+shah_fwd_primers.txt - the guide names and corresponding forward primer sequences
+Shah_metadata_edited.txt - this is an edited copy of the original metadata table
from Shah *et al* with a small transposition error corrected.
+Shah_cut_sites.txt - the locations of the cut sites in the danRer7/Zv9 genome.

The following files refer to Sanger sequencing data from Burger *et al*:

+Burger_Sanger_guides.bed
+Burger_Sanger_metadata.xls
+Burger_wtx_metadata.xls

The following files refer to MiSeq sequencing data from Burger *et al*:

+Burger_MiSeq_guides.bed - this file gives the locations in danRer7 of all guides
used in the Burger *et al* MiSeq sequencing experiment.
+Burger_MiSeq_layout.txt - In Burger *et al*, some PCR primers were used for multiple guides.
This file gives the primer names, the guide names, and the samples the guide occurs in.
+Burger_MiSeq_primer_ranges.bed - this file gives the locations of the PCR primers used in 
the Burger *et al* MiSeq sequencing experiment.

Finally, the following is a copy of the metadata table from Cho *et al*
with some minor reformatting:

+Cho_Table2_reformatted.xls
