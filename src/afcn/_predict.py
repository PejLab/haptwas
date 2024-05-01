"""Predict gene expression from phased genotypes.

By: Genomic Data Modeling Lab
"""

import os
import logging
import numpy as np

from . import model
from . import bedio
from . import vcfio


def run(vcf, par_file, output_fname, filters):

    log2_reference_expression = 0

    if output_fname is None:
        output_fname = "predict"

    output_file = f"{output_fname}.bed"

    logging.info("Begin predictions")
    # open all files
    with (bedio.open_param(par_file, "r") as fpars,
          vcfio.read_vcf(vcf) as fvcf,
          bedio.open_predict(output_file, "w") as fout):

        # write meta data to output_file
        fout.meta["vcf"] = vcf
        fout.meta["parameter_file"] = par_file

        fout.write_meta_data(fvcf.samples)

        # Perform gene expression predictions

        for gene_id, variants in fpars.group_by("gene_id"):

            log2_afc = np.full(len(variants), np.nan)

            haplotypes = [np.full((fvcf.n_samples, len(variants)), np.nan),
                          np.full((fvcf.n_samples, len(variants)), np.nan)]

            # get sample genotypes of each gene associated variant
            for i, v in enumerate(variants):

                # recall that some vcf files have multiple records for a genetic locus.
                # To handle such cases get_genotypes returns a list of
                # records corresponding to that locus, each element of the list
                # is a dictionary of genotype records for each sample.  Most of the
                # time there will be a single record retrieved.
                sample_genotype_records = fvcf.get_genotypes(v[fpars.idx("chrom")],
                                                v[fpars.idx("variant_pos")],
                                                filter_vals=filters)

                # loop over records found for a specific locus.  Break the loop
                # once the parameter file alt allele is a member of the vcf alt alleles

                logging_info = None
                for samp_rec in sample_genotype_records:

                    # what to do with missing data?
                    # alt alleles must match
                    if samp_rec["status"] != 0:
    
                        logging_info = (f"contig:{v[fpars.idx('chrom')]}"
                                        f"\tbed_pos:{v[fpars.idx('variant_pos')]}"
                                        f"\tstatus:{samp_rec['status']}"
                                        f"\tmsg:{samp_rec['msg']}")
                        continue

                    if (alt_allele := v[fpars.idx("alt")]) in samp_rec["alts"]: 
                        logging_info = None
                        break
                    else:
                        logging_info = (f"contig:{v[fpars.idx('chrom')]}"
                                        "\tbed_pos:"
                                        f"{v[fpars.idx('variant_pos')]}"
                                        "\tstatus:-1"
                                        "\tmsg:Alt allele in bed file"
                                        " not a member of alt alleles in VCF.")
                
                # if the loggin_info is not None, than the variant record 
                # retrieval was not successful.  Continue the loop to the next
                # variant
                if logging_info is not None:
                    logging.info(logging_info)
                    continue


                if not samp_rec["phased"]:

                    logging.info(f"contig:{v[fpars.idx('chrom')]}"
                                 f"\tbed_pos:{v[fpars.idx('variant_pos')]}"
                                 "\tstatus:-2"
                                 "\tmsg:Not phased")

                    continue

                log2_afc[i] = v[fpars.idx("log2_afc")]
                for hap_num in range(2):

                    for n, allele_idx in enumerate(samp_rec["genotypes"][hap_num,:]):

                        if allele_idx == 0:
                            haplotypes[hap_num][n, i] = 0
                            continue

                        if allele_idx == samp_rec["genotype_map"][alt_allele]:
                            # alt allele in parameter bed file should be 1

                            haplotypes[hap_num][n, i] = 1


            gene_expr = model.predict(haplotypes[0],
                                      haplotypes[1],
                                      log2_reference_expression,
                                      log2_afc)

            # not all rec values from this iteration of
            # record should
            # have identical genomic coordinates.

            fout.write_line_record(v[fpars.idx("chrom")],
                                   v[fpars.idx("gene_start")],
                                   v[fpars.idx("gene_end")],
                                   gene_id,
                                   gene_expr)

    logging.info("End predictions")
