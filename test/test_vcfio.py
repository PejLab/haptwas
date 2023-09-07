"""Tests of for vcf parsing"""

import os
import unittest
import tempfile

import numpy as np
import pysam

from . import simulate_vcf as simvcf

from afcn import vcfio



def setUpModule():
    global vcf_data 
    global vcf_dir 
    global vcf_name

    vcf_dir = tempfile.TemporaryDirectory()
    vcf_data = simvcf.VCFDataSet()

    tmp_vcf_name = os.path.join(vcf_dir.name, "tmp.vcf")

    with open(tmp_vcf_name, "w") as fid:
        fid.write(str(vcf_data))

    vcf_name = pysam.tabix_index(tmp_vcf_name,
                                 preset="vcf")


def tearDownModule():
    vcf_dir.cleanup()


class TestVCF(unittest.TestCase):
    def setUp(self):
        self.vcf = vcfio.ParseGenotypes(vcf_name, "r")

    def test_properties(self):
        # equal sample size and equal sample names
        self.assertEqual(self.vcf.n_samples,
                         len(vcf_data.samples))

        for test_val,true_val in zip(self.vcf.samples,
                                     vcf_data.samples):
            self.assertEqual(test_val, true_val)

    def test_genotypes_unphased_data(self):
        """Verify genotype records for unphased data."""
        for record in vcf_data.records:

            if record.is_phased():
                continue

            out = self.vcf.get_biallelic_genotypes(record.chrom,
                                                   record.pos-1)

            if record.is_qc_passed() and record.is_biallelic():
                self.assertIsNotNone(out)
            else:
                self.assertIsNone(out)
                continue

            self.assertFalse(out["phased"])

            self.assertEqual(out["alt_allele"],
                             record.alt)

            true_genotypes = np.sum(record.genotypes, 0)

            self.assertTupleEqual(true_genotypes.shape,
                                  out["genotypes"].shape)

            for data_val, true_val in zip(out["genotypes"],
                                          true_genotypes):
                if np.isnan(true_val):
                    self.assertTrue(np.isnan(data_val))
                else:
                    self.assertEqual(data_val, true_val)

    def test_genotype_return_none(self):
        """Verify that None for multiallelic and filtered variants.
        """
        for record in vcf_data.records:

            out = self.vcf.get_biallelic_genotypes(record.chrom,
                                                   record.pos-1)

            if not record.is_biallelic() or record.filter != "PASS":
                self.assertIsNone(out)
            else:
                self.assertIsNotNone(out)


    def test_genotype_phased_data(self):
        """Verify genotype records for phased data."""

        for record in vcf_data.records:
            if not record.is_phased():
                continue

            out = self.vcf.get_biallelic_genotypes(record.chrom,
                                                   record.pos-1)

            if record.is_qc_passed() and record.is_biallelic():
                self.assertIsNotNone(out)
            else:
                self.assertIsNone(out)
                continue

            self.assertTrue(out["phased"])

            self.assertEqual(out["alt_allele"],
                             record.alt)

            true_genotypes = record.genotypes

            self.assertTupleEqual(out["genotypes"].shape, (2, len(self.vcf.samples)))

            for i in range(len(record.samples)):
                if np.isnan(true_genotypes[0, i]):
                    self.assertTrue(np.isnan(out["genotypes"][0,i]))
                else:
                    self.assertEqual(true_genotypes[0, i], out["genotypes"][0, i])

                if np.isnan(true_genotypes[1, i]):
                    self.assertTrue(np.isnan(out["genotypes"][1,i]))
                else:
                    self.assertEqual(true_genotypes[1, i], out["genotypes"][1, i])

    def test_no_variants_foundo(self):
        """Verify genotype records for phased data."""
        # note that in the simulated vcf there is no variant
        # at position 0
    
        for record in vcf_data.records:
            break

        out = self.vcf.get_biallelic_genotypes(record.chrom, 1)
        self.assertIsNone(out)
        

if __name__ == "__main__":
    unittest.addModuleCleanup(tearDownModule)
    unittest.main()
    unittest.doModuleCleanup()
