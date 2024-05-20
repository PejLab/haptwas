"""Test gene / variant map

By: GDML
"""

import os
import shutil
import io
import tempfile
import time

import unittest

from afcn import bedio


class TestFactoryFunction(unittest.TestCase):
    def test_wrong_file_extensions(self):
        with self.assertRaises(NotImplementedError):
            bedio.read_gene_variant_map("test.vcf")

        with self.assertRaises(NotImplementedError):
            bedio.read_gene_variant_map("test.bgz")




if __name__ == "__main__":
    unittest.main()

