# MOCHI

## Required Package
MOCHI requires:
* Python (tested 3.6.8)
* numpy (tested 1.13.1)
* pandas (tested 0.20.3)


## Usage
The parameters are as followed:

Cd Code

python entrance.py -c gm12878 -m 4

* "-c","--cell",default = 'gm12878' :Controls the data it runs
* "-p", "--parent_threshold",default = 0.45 :Controls the threshold of the overlapping TF
* "-t", "--total",default = False :Whether to run the code all over again. By default it uses the cached data to accelerate.
* "-m", "--motif",default = 4 :4-node motif or 3-node motif



## Input Example

An input example for cell line gm12878 is provided at the directory ./Data/gm12878/used/

Required files are:

**hic_gm12878_KR_10000_intra_subset_chr*_oe_1_genes.txt**

Intra chromosomal matrix with KR and oe normalization. The data is filtered to kept only the region with expressed genes:

The file looks like:

chr1:890000 chr1:890000 3459.2085
chr1:890000 chr1:900000 1114.036
chr1:900000 chr1:900000 3740.9792

...

**hic_gm12878_KR_10000_inter_top_1_genes.txt**

Inter chromosomal matrix with KR normalization. The data is filtered to kept only the region with expressed genes:

The file looks like:

chr2:114380000 chr22:51220000 1288.5735
chr1:109650000 chr22:30160000 105.81676
chr7:26240000 chr15:40850000 87.49172

...

**grn_gm12878_txStart.txt**

Gene regulatory network with gene names.

The file looks like:

Source	Target
CBX5	FAM53B
ZBTB49	CUTC
ETV3	BYSL

...

**gene_chrom_bin_num_hg19_gm12878.txt**

A file to map gene name to its corresponding bin coordinates.

The file should at least have the following columns:

Chrom	Start	end	Gene_name(TF)	Bin_10kb

chr12	9220303	9268558	A2M	chr12:9260000

...

The column "Gene_name" is named as "TF in the example input.

**Extra file for analysis** 

To repeat some of the analysis we did in the manuscript, annotations like PPI network, compartments, loop, repli-seq... should be provided.

We provide the corresponding files in the directory  ./Statistics

