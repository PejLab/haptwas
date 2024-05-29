""""Tests for the prediction bed file IO.

By: GDML
"""
import os
from datetime import datetime
import io
from unittest import TestCase, main

from afcn import bedio, utils


class GeneRecord:
    _field_delimiter = "\t"
    _phase_delimiter = "|"
    _new_line = "\n"

    def __init__(self,chrom, start, end, gene_id, h1, h2):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.gene_id = gene_id
        self.h1 = h1
        self.h2 = h2
 
    def __repr__(self):
        s = self._field_delimiter.join([self.chrom,
                                        str(self.start),
                                        str(self.end),
                                        self.gene_id])

        for h1, h2 in zip(self.h1, self.h2):
            s += self._phase_delimiter.join([f"{self._field_delimiter}{h1}",
                                             str(h2)])

        return s


class DataRecords:
    def __init__(self):
        self._records = []

    def __iter__(self):
        for rec in self._records:
            yield rec

    def append(self, chrom, start, end, gene_id, h1, h2):
        self._records.append(GeneRecord(chrom, start, end, gene_id, h1, h2))

    def get_record_by_idx(self, idx):
        return self._records[idx]


class TestWritePredictionBed(TestCase):

    encoding = "utf-8"
    _py_version=".".join([str(os.sys.version_info.major),
                          str(os.sys.version_info.minor),
                          str(os.sys.version_info.micro)])

    if "USER" in os.environ:
        user = os.environ["USER"]
    elif "LOGIN" in os.environ:
        user = os.environ["LOGIN"]
    else:
        user = "dont_know"

    meta_prefix = "##"
    header_prefix = "#"

    # meta data
    date = datetime.today()
    meta_data = [f"{meta_prefix}afcn_version={utils.get_version()}",
                 f"{meta_prefix}date={date.year}-{date.month:02d}-{date.day:02d}",
                 f"{meta_prefix}python_version={_py_version}",
                 f"{meta_prefix}user={user}"]

    sample_names = ["sample_1", "sample_2"]
    header = "\t".join([f"{header_prefix}chrom",
                        "start",
                        "end",
                        "gene_id",
                        *sample_names])

    h1_preds = [[23.2, 1],
                [2,21]]
    h2_preds = [[2, 21],
                [3, 1.543]]

    i = 1

    records = DataRecords()

    for h1, h2 in zip(h1_preds, h2_preds):
        records.append("chr1",
                        1000 * i, 1304 * i,
                        f"ENSG00FAKEGENE{i}",
                        h1, h2)

        
    def setUp(self):
        """Make a list of target prediction file lines"""
        self.b_id = io.BytesIO()

        self.target_predict_bed_lines = []

        for tline in self.meta_data:
            self.target_predict_bed_lines.append(tline)

        self.target_predict_bed_lines.append(self.header) 

        for record in self.records:
            self.target_predict_bed_lines.append(str(record))


    def tearDown(self):
        if not self.b_id.closed:
            self.b_id.close()


    def test_init(self):
        fid = bedio.WritePredictionBed(self.b_id)

        self.assertEqual(fid._fid, self.b_id)

        self.assertFalse(fid._meta_and_header_written)

        self.b_id.close()
        self.assertTrue(fid.closed)

    def test_context_management(self):
        fid = bedio.WritePredictionBed(self.b_id)

        self.assertFalse(fid.closed)
        self.assertFalse(self.b_id.closed)

        self.assertEqual(fid, fid.__enter__())
        fid.__exit__()

        self.assertTrue(fid.closed)
        self.assertTrue(self.b_id.closed)
        
    def test_write_meta_data(self):
        with bedio.WritePredictionBed(self.b_id) as fid:

            fid.write_meta_data(self.sample_names)
        
            fid._fid.seek(0)

            for correct_line in self.meta_data:
                rline = fid._fid.readline().decode(self.encoding).strip()
                self.assertEqual(rline,
                                 correct_line)

            self.assertTrue(fid._meta_and_header_written)
        
    def test_exeception(self):
        with (bedio.WritePredictionBed(self.b_id) as fid,
              self.assertRaises(ValueError)):
            fid.write_line_record("1", 10, 100, "SIM_GENE_1",
                                  3, 1.5)

    def test_header(self):

        with bedio.WritePredictionBed(self.b_id) as fid:
            fid.write_meta_data(self.sample_names)

            fid._fid.seek(0)
            for tline in fid._fid:
                pass

        self.assertEqual(tline.decode(encoding=self.encoding).strip(),
                         self.header)

    def test_records(self):
        
        with bedio.WritePredictionBed(self.b_id) as fid:
            fid.write_meta_data(self.sample_names)

            for r in self.records:
                
                fid.write_line_record(r.chrom,
                                      r.start,
                                      r.end,
                                      r.gene_id,
                                      r.h1,
                                      r.h2)

            fid._fid.seek(0)

            i = 0
            for tline in fid._fid:

                tline = tline.decode(encoding=self.encoding).strip()

                if tline.startswith(self.header_prefix):
                    continue

        
                self.assertEqual(tline, str(self.records.get_record_by_idx(i)))

                i += 1






if __name__ == "__main__":
    main()
