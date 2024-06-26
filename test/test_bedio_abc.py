""""Test for bedio abstract base classes

By: Robert Vogel
    Genomic Data Modeling Lab
"""
import os
from datetime import datetime
import tempfile
import io
from collections import OrderedDict
from unittest import TestCase, main

from afcn import bedio, utils


class TestBedABC(TestCase):
    def setUp(self):
        _file_suffix = ".bed"
        self.encoding = "utf-8"

        _meta_prefix = "##"
        _meta_data_key_val_delimiter = "="

        _header_prefix = "#"

        _field_delimiter = "\t"
        _new_line = "\n"

        _date = datetime.today()
        _date = f"{_date.year}-{_date.month:02d}-{_date.day:02d}"

        _py_version=".".join([str(os.sys.version_info.major),
                              str(os.sys.version_info.minor),
                              str(os.sys.version_info.micro)])

        if "USER" in os.environ:
            _user = os.environ["USER"]
        elif "LOGIN" in os.environ:
            _user = os.environ["LOGIN"]
        else:
            _user = "not_available"

        self.meta = OrderedDict(afcn_version=utils.get_version(),
                                date = _date,
                                python_version = _py_version,
                                user = _user)

    def test_bed_spec(self):
        b = bedio.BedABC()
        self.assertEqual(b._file_suffix , ".bed")

        self.assertEqual(b._encoding , "utf-8")

        self.assertEqual(b._meta_prefix , "##")
        self.assertEqual(b._meta_data_key_val_delimiter , "=")

        self.assertEqual(b._header_prefix , "#")

        self.assertEqual(b._field_delimiter , "\t")
        self.assertEqual(b._new_line , "\n")


    def test_initialize(self):
        b = bedio.BedABC()

        self.assertIsNone(b._fid)
        self.assertIsNone(b.header)

        meta_data_keys = ("afcn_version", "date", "python_version", "user")

        for key in meta_data_keys:
            self.assertIn(key, self.meta)

        self.assertEqual(len(self.meta), len(meta_data_keys))

        self.assertIsInstance(b.meta, OrderedDict)




class TestParseBedABC(TestCase):
    def setUp(self):
        self.encoding="utf-8"

        # make in memory file to parse
        b = bedio.BedABC()
        self._fid = io.BytesIO()

        for key, val in b.meta.items():
            self._fid.write((f"{b._meta_prefix}{key}{b._meta_data_key_val_delimiter}"
                            f"{val}{b._new_line}").encode(self.encoding))

        self.header = ["chrom","start", "end", "name"]
        self.req_header_fields = dict(chrom=str, start=int, end=int, name=str)

        self._fid.write("".join([b._header_prefix,
                    b._field_delimiter.join(self.header)]).encode(self.encoding))


        self.records = [["chrm1", "1000", "1100", "ENS_TEST_002"],
                         ["chrm1", "1200", "1300", "ENS_TEST_002"],
                         ["chrm1", "1600", "1900", "ENS_TEST_003"]]

        for rec in self.records:
            self._fid.write("".join([b._new_line,
                             b._field_delimiter.join(rec)]).encode(self.encoding))


        self._fid.seek(0)

        self.parser = bedio.ParseBedABC(self._fid)


    def tearDown(self):
        if not self._fid.closed:
            self._fid.close()

    def test_initialize(self):
        self.assertEqual(self._fid, self.parser._fid)

        self.assertIsInstance(self.parser._colname_to_idx, dict)

    def test_file_operation_property(self):
        self.assertFalse(self.parser.closed)
        self.parser.close()
        self.assertTrue(self.parser.closed)
        self.assertTrue(self._fid.closed)

    def test_context_management(self):
        with bedio.ParseBedABC(self._fid) as fid:
            self.assertFalse(fid.closed)

        self.assertTrue(fid.closed)

    def test_record_parser(self):
        with self.assertRaises(NotImplementedError):
            self.parser._record_parser()

    def test_colname_error(self):
        with self.assertRaises(KeyError):
            self.parser.idx("cat")

    def test_colname(self):
        self.parser._req_header_fields = self.req_header_fields
        self.parser._initialize()

        for i, colname in enumerate(self.header):
            self.assertEqual(self.parser.idx(colname), i)

        self.assertEqual(len(self.parser._colname_to_idx),
                         len(self.header))

    def test_initialize(self):
        self.parser._req_header_fields = self.req_header_fields
        self.parser._initialize()

        for i, colname in enumerate(self.header):
            self.assertEqual(colname, self.parser.header[i])

    def test_empty_file(self):
        with (self.assertRaises(UnboundLocalError),
            bedio.ParseBedABC(io.BytesIO()) as parser):
            parser._initialize()


    def _record_parser(self):
        for rec in self.parser:
            output = (rec.strip()
                        .split(self.parser._field_delimiter))

            i = 0

            for _, field_type in self.parser._req_header_fields.items():
                output[i] = field_type(output[i])
                i += 1

            yield output

    def test_group_by(self):
        self.parser._record_parser = self._record_parser
        self.parser._req_header_fields = self.req_header_fields
        self.parser._initialize()

        i = 0
        for gene_name, g in self.parser.group_by("name"):
            print(gene_name, g)
            self.assertEqual(gene_name, self.records[i][-1])

            for true_record, test_record in zip(self.records[i:len(g)], g):

                j = 0
                for rec, val in zip(true_record, test_record):
                    field_type = self.req_header_fields[self.header[j]]
                
                    rec = field_type(rec)
                    self.assertEqual(rec, val)
                    j += 1

                i += 1


if __name__ == "__main__":
    unittest.addModuleCleanup(tearDownModule)
    main()
    unittest.doModuleCleanup()
