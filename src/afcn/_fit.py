
import os
import logging
import numpy as np


def run(vcf, gene_expr, gene_and_variants, out_dir):
    logging.info("Start: effect size inference")
    logging.info(f"vcf file:{vcf}")
    logging.info(f"expression file:{gene_expr}")

    print(f"Successfully found {__file__}.")

    logging.info(f"output written to:{out_dir}/predictions.bed")
    logging.info("Finished")

    raise NotImplementedError
