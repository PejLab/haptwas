"""VCF IO tools.

By: GDML
"""
import os
import numpy as np
from pysam import VariantFile


# TODO automatic generation of tabix index from bgzip vcf
def read_vcf(vcf):
    """Open file object """
    tbx_file = f"{vcf}.tbi"
    if not os.path.exists(tbx_file):
        raise ValueError("Need tabix index for vcf file.")

    if os.path.exists(tbx_file):
        return ParseGenotypes(vcf, "r", index_filename=tbx_file)

    raise RuntimeWarning("No tabix file detected.")


# TODO confirm indexing rules.
class ParseGenotypes(VariantFile):
    @property
    def samples(self):
        return list(self.header.samples)

    @property
    def n_samples(self):
        return len(self.header.samples)

    @property
    def contigs(self):
        return list(self.header.contigs)

    def len_contig(self, contig):
        return self.header.contigs.get(contig).length

    def get_contig_variants(self, contig):
        positions = dict()

        for record in self.fetch(contig):
            positions[record.id] = record.pos

        return positions

    def get_biallelic_genotypes(self, contig, pos, filter_val="PASS"):
        """Get array(s) of genotypes of a biallelic locus.

        Extract the alt allele count, unphased {0,1,2,np.nan} or phased
        {0,1,np.nan} for all samples of at a specified locus.  Only
        return data for biallelic variants.  For a genotypes of a variant
        to be considered phased, all samples must be phased.

        Args:
            contig: (str)
                the contig in which the variant is located.
            pos: (int)
                assume 0 based indexing

        Returns:
            dict: 
                {
                    phased: (bool),
                    genotypes: ((n_sample,) np.ndarray) if not phased,
                        otherwise ((2, n_sample) np.ndarray)
                    alt_allele: (character)
                    status: (int)
                    msg: (str)
                }

        where `status` and `msg` reports whether a variant that satisfies
        our requirements is found, if not it reports the incompatibility.

        | status code    |  msg                           |
        | -------------- | ------------------------------ |
        | 0              | Success                        |
        | 1              | More than one variant found.   |
        | 2              | No variant found.              |
        | 3              | No alt allele found.           |
        | 4              | More than 1 alt allele found.  |
        | 5              | Filter mismatch                |
        """

        genotypes = np.zeros(shape=(2, self.n_samples))
        phased = True
        
        # Set i to None to detect when no variant is found.  When no
        # variant is found the for loop below does not assign a value
        # to i.
        i = None

        for i, variant in enumerate(self.fetch(contig, pos, pos+1)):
            if i > 0:
                return dict(status= 1,
                            msg="More than one variant found")

        # check whether a variant was found
        if i is None:
            return dict(status=2,
                        msg="No variant found")

        # return none if variant isn't an SNV, and doesn't pass
        # filter, capitilization matters:
        #   * deletion, e.g. Ref field = "."
        #   * multiallelic variant
        #   * more than one filter value
        #   * specified filter satsified
        if variant.alts is None:
            return dict(status=3,
                        msg="No alt allele")

        if len(variant.alts) > 1:
            return dict(status=4,
                        msg="More than 1 alt allele found")

        if (len(filter_vals := variant.filter.keys()) != 1
            or filter_val not in filter_vals):

            return dict(status=5,
                        msg=("VCF filter value(s) does exclusively equal"
                             f" the required input of {filter_val}"))


        # Note that pysam returns None for missing alleles in a sample
        var_encoding = {variant.ref:0,
                        variant.alts[0]:1,
                        None:np.nan}

        for n, samp_geno in enumerate(variant.samples.itervalues()):

            genotypes[0, n] = var_encoding[samp_geno.alleles[0]]
            genotypes[1, n] = var_encoding[samp_geno.alleles[1]]

            # if geno.phased is False for any one sample, then by
            # this if statement geno.phased will be False at the
            # termination of for loop
            if phased:
                phased = samp_geno.phased

        if not phased:
            genotypes = np.sum(genotypes, 0)

        return {"phased":phased,
                "genotypes":genotypes,
                "alt_allele": variant.alts[0],
                "status": 0,
                "msg": "success"}
