"""Test gene / variant map

By: Robert Vogel
    Genomic Data Modeling Lab
"""

import os
import shutil
import io
import tempfile
import time

import unittest

from collections import OrderedDict
from afcn import bedio


class TestFactoryFunction(unittest.TestCase):
    def test_wrong_file_extensions(self):
        with self.assertRaises(NotImplementedError):
            bedio.read_eqtl_map("test.vcf")

        with self.assertRaises(NotImplementedError):
            bedio.read_eqtl_map("test.bgz")


class TestParseEqtlMap(unittest.TestCase):
    @staticmethod
    def make_meta_record(key, val, newline=False):
        s = (f"{bedio.BedABC._meta_prefix}{key}"
             f"{bedio.BedABC._meta_data_key_val_delimiter}"
             f"{val}")

        if newline:
            return f"{bedio.BedABC._new_line}{s}".encode(bedio.BedABC._encoding)

        return s.encode(bedio.BedABC._encoding)

    def make_header_str(self, newline=True):
        
        s = "".join([bedio.BedABC._header_prefix,
                    bedio.BedABC._field_delimiter.join(self.req_header.keys())])

        if newline:
            return f"{bedio.BedABC._new_line}{s}".encode(bedio.BedABC._encoding)

        return s.encode(bedio.BedABC._encoding)

    def make_data_record(self, rec):
        # chrom, a string is the first item
        s = f"{bedio.BedABC._new_line}{rec[0]}"

        for field_rec in rec[1:]:
            s+= f"{bedio.BedABC._field_delimiter}{str(field_rec)}"

        return s.encode(bedio.BedABC._encoding)

    def setUp(self):
        self.meta = {"user":"me"}

        self.req_header = OrderedDict(
                    chrom = str,
                    qtl_start = int,
                    qtl_end = int,
                    qtl_id = str,
                    gene_start = int,
                    gene_end = int,
                    gene_id = str,
                    variant_pos = int,
                    variant_id = str,
                    ref = str,
                    alt = str)

        self.header = list(self.req_header.keys())

        self.records = [
                ["chrm1",10, 1000,"qtl_id_1",
                 20, 1000, "ENSG_1", 10, "var_1_id", "A","T"],
                ["chrm1",100, 10000,"qtl_id_2",
                 200, 10000, "ENSG_1", 100, "var_2_id", "G","C"],
                ["chrm1",100, 10000,"qtl_id_3",
                 200, 10000, "ENSG_2", 100, "var_3_id", "G","C"]
                ]

        self.gene_idx = 6
        
        self._fid = io.BytesIO()

        self._fid.write(self.make_meta_record("user", self.meta["user"]))
        self._fid.write(self.make_header_str())

        for rec in self.records:
            self._fid.write(self.make_data_record(rec))

        self._fid.seek(0)

        self.parser = bedio.ParseEqtlMap(self._fid)

    def tearDown(self):
        self._fid.close()

    def test_fid(self):
        self.assertEqual(self._fid, self.parser._fid)

    def test_group_by(self):

        i = 0
        for gene_name, g in self.parser.group_by("gene_id"):

            self.assertEqual(gene_name, self.records[i][self.gene_idx])

            for true_record, test_record in zip(self.records[i:len(g)], g):

                j = 0
                for rec, val in zip(true_record, test_record):
                    field_type = self.req_header[self.header[j]]
                
                    rec = field_type(rec)
                    self.assertEqual(rec, val)
                    j += 1

                i += 1


if __name__ == "__main__":
    unittest.main()

