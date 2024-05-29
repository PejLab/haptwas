### SPECIFICATION: GENE EXPRESSION BED

BED 4 file with meta data and header as follows:

* line [0, n): meta data prepended with ##
* line n: header prepended with #
* line [n+1, n+1+N_samples): records

Contents:

* meta data
    - ##afcn_version=(version str)
    - ##vcf_file=(str)
    - ##parameter_file=(str)
* BED fields labeled: 
    - chrom
    - start : gene start position
    - end : gene end position
    - name : ensembl gene id
* custom fields:
    - sample_id: (float) predicted gene expression per haplotype, e.g.
        322|32
