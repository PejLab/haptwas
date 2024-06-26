"""Fit afc model to data

By: Genomic Data Modeling Lab
"""

import os
import logging
import numpy as np

from . import model
from . import bedio
from . import vcfio


def run(vcf, gene_expr, eqtls, output_path_and_prefix):

    output_fname = f"{output_path_and_prefix}{os.path.extsep}bed"

    logging.info(f"output_file:{output_fname}")
    logging.info("Begin effect size inference")

    with (vcfio.read_vcf(vcf) as fin_vcf,
          bedio.read_gene_expression(gene_expr) as fin_expr,
          bedio.read_gene_variant_map(eqtls) as fin_eqtls,
          bedio.open_param(output_fname, "w") as fout_param):
        pass

     
    logging.info("Finished")

    raise NotImplementedError
