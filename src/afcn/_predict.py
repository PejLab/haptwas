"""Predict gene expression from phased genotypes.

By: Genomic Data Modeling Lab
"""

import os
import numpy as np

from . import model
from . import bedio
from . import vcfio


def run(vcf, par_file, output_dir, reference_expression=0):

    log2_reference_expression = 0

    output_file = os.path.join(output_dir, "predictions.bed")

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

            hap_one = np.full((fvcf.n_samples, len(variants)), np.nan)
            hap_two = np.full((fvcf.n_samples, len(variants)), np.nan)

            # get sample genotypes of each gene associated variant
            for i, v in enumerate(variants):
                tmp = fvcf.get_biallelic_genotypes(v[fpars.idx("chrom")],
                                                   v[fpars.idx("variant_pos")])

                # what to do with missing data?
                # alt alleles must match
                if (tmp is None
                    or tmp["alt_allele"] != v[fpars.idx("alt")]):
                    continue
                
                # TODO should I raise exception or just not include
                # variant in analysis?
                if not tmp["phased"]:
                    continue
                #raise ValueError(("Data are not phased, "
                #                    "as is required."))

                log2_afc[i] = v[fpars.idx("log2_afc")]

                hap_one[:, i] = tmp["genotypes"][0, :]
                hap_two[:, i] = tmp["genotypes"][1, :]

            gene_expr = model.predict(hap_one,
                                      hap_two,
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
