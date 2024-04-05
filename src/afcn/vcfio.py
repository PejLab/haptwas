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

    def get_genotypes(self, contig, pos, filter_vals=["PASS"]):
        """Get array(s) of genotypes.

        Extract the genotype{0,1,2,...,np.nan} for all samples
        of at a specified locus. For a genotypes of a variant
        to be considered phased, all samples must be phased.

        Args:
            contig: (str)
                the contig in which the variant is located.
            pos: (int)
                assume 0 based indexing
            filter_vals:(list)
                list of strings specifing acceptable
                variant filter value(s)

        Returns:
            dict: 
                {
                    phased: (bool),
                    genotypes: ((2, n_sample) np.ndarray)
                    var_encodings: (dict)
                        key value pairs of allele (A,T,C,G)
                        and numeric value.  Missing genotypes
                        are None
                    status: (int)
                    msg: (str)
                }

        where `status` and `msg` reports whether a variant that
        satisfies our requirements is found, if not it reports
        the incompatibility.

        | status code    |  msg                           |
        | -------------- | ------------------------------ |
        | 0              | Success                        |
        | 1              | More than one variant found.   |
        | 2              | No variant found.              |
        | 3              | No alt allele found.           |
        | 4              | Filter mismatch                |
        """

        if not isinstance(filter_vals, list):
            raise TypeError("Input filter value must be a Python list")

        genotypes = np.full((2, self.n_samples), np.nan)
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
        #   * more than one filter value
        #   * specified filter satsified
        if variant.alts is None:
            return dict(status=3,
                        msg="No alt allele")

        # the lower operation ensures that pattern matching
        # is not case sensitive
        satisfy_filter_criterion = False
        variant_filter_vals = []
        for w in variant.filter.keys():
            variant_filter_vals.append(w.lower())

        if len(variant_filter_vals) == 0 :
            variant_filter_vals.append("missing")

        for fval in filter_vals:
            if fval.lower() in variant_filter_vals:
                satisfy_filter_criterion = True
                break

        if not satisfy_filter_criterion:
            return dict(status=4,
                        msg=("VCF filter value(s) does not match"
                             " required input of"
                             f" ({','.join(filter_vals)})"))


        # Note that pysam returns None for missing alleles in
        # a sample

        # according to vcf4.4 specification
        genotype_map = {variant.ref:0,
                        None:np.nan}


        for alt_allele in variant.alts:
            genotype_map[alt_allele] = None


        for n, samp_geno in enumerate(variant.samples.itervalues()):

            # sets the alternative allele encodings
            for allele, idx in zip(samp_geno.alleles,
                                   samp_geno.allele_indices):


                # default entry is np.nan
                if allele is None and idx is None:
                    continue
                elif genotype_map[allele] is None:
                    genotype_map[allele] = idx
                elif genotype_map[allele] != idx:
                    raise ValueError("Inconsistent indexing of"
                                     " alleles across samples")


            genotypes[0, n] = genotype_map[samp_geno.alleles[0]]
            genotypes[1, n] = genotype_map[samp_geno.alleles[1]]

            # if geno.phased is False for any one sample, then by
            # this if statement geno.phased will be False at the
            # termination of for loop
            if phased:
                phased = samp_geno.phased

        return {"phased":phased,
                "genotypes":genotypes,
                "genotype_map":genotype_map,
                "alts":variant.alts,
                "status": 0,
                "msg": "success"}
