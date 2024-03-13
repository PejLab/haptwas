"""Predict gene expression from phased genotypes.

By: Genomic Data Modeling Lab
"""

import os
import logging
import numpy as np

from . import model
from . import bedio
from . import vcfio


def run(vcf, par_file, output_fname, reference_expression=0):

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
                tmp = fvcf.get_genotypes(v[fpars.idx("chrom")],
                                    v[fpars.idx("variant_pos")])

                # what to do with missing data?
                # alt alleles must match
                if tmp["status"] != 0:

                    logging.info(f"contig:{v[fpars.idx('chrom')]}"
                                 f"\tbed_pos:{v[fpars.idx('variant_pos')]}"
                                 f"\tstatus:{tmp['status']}"
                                 f"\tmsg:{tmp['msg']}")
                    continue

                if (alt_allele := v[fpars.idx("alt")]) not in tmp["alts"]: 

                    logging.info(f"contig:{v[fpars.idx('chrom')]}"
                                 "\tbed_pos:"
                                 f"{v[fpars.idx('variant_pos')]}"
                                 "\tstatus:-1"
                                 "\tmsg:Alt allele in bed file"
                                 " not a member of alt alleles in VCF.")
                    continue
                
                # TODO should I raise exception or just not include
                # variant in analysis?
                if not tmp["phased"]:

                    logging.info(f"contig:{v[fpars.idx('chrom')]}"
                                 f"\tbed_pos:{v[fpars.idx('variant_pos')]}"
                                 "\tstatus:-2"
                                 "\tmsg:Not phased")

                    continue

                log2_afc[i] = v[fpars.idx("log2_afc")]

                for hap_num in range(2):

                    for n, allele_idx in enumerate(tmp["genotypes"][hap_num,:]):

                        if allele_idx == 0:
                            haplotypes[hap_num][n, i] = 0
                            continue

                        if allele_idx == tmp["genotype_map"][alt_allele]:
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
