# HotDriver
HotDriver is a python program that is designed to employ the Poisson statistical model to predict hotspot mutations in single amino acid level given a mutation dataset.

Prerequisites for running HotDriver:
1. python 2.7.2 or above;
2. additional python package: scipy

In /tools, bedtools (a binary version) is provided for the userÕs convenience, the user could also choose to use another version of bedtools, but need to change the path of bedtools in line 264 of HotDriver.py.

In /example_data, multiple examples of data are provided.
1) CpG island bed file, in which the CpG regions in human genome (hg19 UCSC) are defined;
Format: 
chromosome	start	end

2) gene_length.tsv, in which the genes and their corresponding amino acid lengths are defined;
Format: 
gene	AA_length

3) gene_list.tsv, in which a list of candidate genes are provided;
Format:
gene

4) mutation_data.tsv, in which the mutation data that can be used to predict hotspot mutations are provided.
Format:
Sample_ID	Gene_ID	Chromosome	Start	End	Ref_allele	Alt_allele	Protein_variant	Mutation_subtype
	
The user could choose different data sources, but the format of each individual file has to be exactly the same as that of the example files to make sure HotDriver work properly.

Please contact tchen1@mdanderson.org or kchen3@mdanderson.org if the user may have any questions.
