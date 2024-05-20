"""IO tools for afcn bed files.

By: GDML
"""

import os
import io
import itertools
from collections import OrderedDict
import gzip
from datetime import datetime

import numpy as np

from . import utils

# TODO how to manage specification version, e.g. code may 
# change but spec does not.
# TODO where to write specification so that it is callable 
# to __main__.py and this code

BED_SUFFIX = ".bed"

def open_param(filename, mode):
    """Open either gzipped or normal text file for reading.

    Args:
        filename: (str)
            either .bed or .bed.gz file
        mode: (char)
            "r" for read, "w" for write

    Returns:
        returns instance of class that mimics / uses io.FileIO
        services
    """

    if (not filename.endswith(".bed.gz") and 
        not filename.endswith(".bed")):

        raise ValueError("Not a bed file.")

    if not os.path.exists(filename):
        raise FileNotFoundError(filename)

    if mode == "r" and filename.endswith(".bed"):
        return ParseParamBed(io.FileIO(filename, mode))

    elif mode == "r" and filename.endswith(".bed.gz"):
        return ParseParamBed(gzip.GzipFile(filename))

    if mode == "w" and filename.endswith(".bed"):
        return WriteParamBed(io.FileIO(filename, mode))

    raise NotImplementedError


class BedSpec:
    _encoding = "utf-8"

    _meta_prefix = "##"
    _meta_data_key_val_delimiter = "="

    _header_prefix = "#"

    _field_delimiter = "\t"
    _new_line = "\n"

    def __init__(self, *_):
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
        self.header = None
        self._colname_to_idx  = dict()

    def __iter__(self):
        return self

    def __next__(self):
        return self._fid.__next__().decode(self._encoding)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        if not self.closed:
            self.close()

    def tell(self):
        return self._fid.tell()

    def seek(self, *args):
        return self._fid.seek(*args)

    @property
    def closed(self):
        return self._fid.closed

    def close(self):
        if not self.closed:
            self._fid.flush()
            self._fid.close()


class ExpressionQTLBedSpec(BedSpec):
    _req_header_fields = OrderedDict(chrom = str,
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


class ParamBedABC(ExpressionQTLBedSpec):
    """Enforce parameter bed file specification."""

    def __init__(self):
        super().__init__()

        self._req_header_fields["log2_afc"] = float
        self._req_header_fields["sem"] = float
        self._req_header_fields["p_val"] = float

        for idx, col_name in enumerate(self._req_header_fields.keys()):
            self._colname_to_idx[col_name] = idx


class ParseParamBed(ParamBedABC):
    def __init__(self, fid):

        self._fid = fid

        # recall that the order of class methods being called is:
        #   1. __init__
        #   2. __enter__
        #   3. __exit__
        # If there exists an error in __init__ or __enter__, the __exit__
        # method is not guaranteed to run.  Consequently, I use the try block.
        try:

            super().__init__()
            self._initialize()

        except Exception:

            self.__exit__()
            raise

    def _initialize(self):
        """Load and validate meta data and header as defined in spec."""

        # load meta data 
        for par_line in self:

            if not par_line.startswith(self._meta_prefix):
                break

            # max_split kwarg in the string split method is important in 
            # the event that a user has an '=' in the value field.

            key, val = (par_line
                        .removeprefix(self._meta_prefix)
                        .split(sep=self._meta_data_key_val_delimiter,
                               maxsplit=1))

            self.meta[key.strip()] = val.strip()

        # decompose header and verify it is spec. compliant
        # remember that for empty files par_line is not associated 
        # with any value and throws an UnboundLocalError
        if not par_line.startswith(self._header_prefix):
            raise ValueError("Header required: No file header found")


        try:
            
            par_line = par_line.removeprefix(self._header_prefix)

        except UnboundLocalError as err:
            # the add_note method is available on Python 3.11 +
            #err.add_note(("\nMost likely reason for failure is that "
            #              f"file {self.name} is empty."))
            err.args = (err.args[0] + 
                        f"\nMost likely reason for failure is "
                        "that file is empty.",)
            raise UnboundLocalError(err) from None

        self.header = par_line.strip().split(sep=self._field_delimiter)

        for i, field_name in enumerate(self._req_header_fields.keys()):

            if field_name != self.header[i]:
                raise ValueError(f"Invalid file header")

        self._data_char_number = self.tell()

    def _record_parser(self):

        if self.header is None:
            return None

        for record in self:

            output = record.strip().split(self._field_delimiter)

            i = 0
            for _, field_type in self._req_header_fields.items():
                output[i] = field_type(output[i])
                i += 1

            for out_val in output[i:]:

                if utils.is_int(out_val):
                    out_val = int(out_val)
                elif utils.is_float(out_val):
                    out_val = float(out_val)

                output[i] = out_val

                i += 1

            yield output

    def idx(self, col_name):
        if col_name not in self._colname_to_idx:
            raise KeyError(f"{col_name} is not in header")
        return self._colname_to_idx[col_name]

    def group_by(self, gb_id):
        """Perform group by on gb_id, assumed sorted.

        Requires that gene id are sorted lexicographically.
        """

        if (self.tell() != self._data_char_number):
            self.seek(self._data_char_number)

        # start with ! as it has lowest lexicographic order
        previous_record_gb_id = "!"

        for k, g in itertools.groupby(self._record_parser(), 
                                      key=lambda x: x[self.idx(gb_id)]):

            if k < previous_record_gb_id:
                raise ValueError(f"{gb_id} column is not sorted.")

            previous_record_gb_id = k

            yield k, list(g)


# TODO
class WriteParamBed(ParamBedABC):
    def __init__(self, fid):

        self._fid = fid
        
        try:

            super().__init__()

            self._meta_and_header_written = False

            self.meta["data"] = "Inferred afc model parameters."

        except Exception:

            self.__exit__()
            raise


    def write_meta_and_header_data(self, sample_names):
        # instantiate string
        s = ""

        # construct meta data
        for key, val in self.meta.items():
            s += "".join([self._meta_prefix,
                          key,
                          self._meta_data_key_val_delimiter,
                          val,
                          self._new_line])

        # construct header
        s += self._header_prefix
        s += self._field_delimiter.join(*self._req_header_fields.keys())
        s = s.encode(encoding=self._encoding)

        if self._fid.write(s) != len(s):
            raise RuntimeError("Bytes written not equal to bytes of string.")

        self._meta_and_header_written = True

    def write_line_record(self):
        """Write line of predictions.

        Args:
            chrom: (str) chromosome name
            qtl_start: (int) genomic coordinate of gene beginning
            qtl_end: (int) genomic coordinate of gene end
            qtl_id: (str)
            gene_star: (int)
            gene_end: (int)
            gene_id: (int) Ensembl gene id
            variant_pos: (int) Single nucleotide polymorphism
            variant_id: (str),
            ref: (str), reference allele
            alt: (str), alternative allele
            log2_afc: (float)
            log2_afc_sem: (float)
            log2_afc_p_value: (float)


        Returns:
            None
        """
        pass
#         if not self._meta_and_header_written:
#             raise ValueError("Write meta data before real data")
# 
#         data_str = self._field_delimiter.join([str(w) for w in data])
#         chrom = self._new_line + chrom
# 
#         record_str = self._field_delimiter.join([chrom, 
#                                                 str(start), 
#                                                 str(end), 
#                                                 name,
#                                                 data_str])
#         record_str = record_str.encode(encoding=self._encoding)
# 
#         if self._fid.write(record_str) != len(record_str):
#             raise RuntimeError("Bytes written not equal to bytes of string.")

def read_gene_variant_map(filename):
    if filename.endswith(".bed"):
        return ParseGeneVariant(io.FileIO(filename))

    if filename.endswith(".bed.gz"):
        return ParseGeneVariantMap(gzip.GzipFile(filename))

    raise NotImplementedError


# TODO

class ParseGeneVariantMap(ExpressionQTLBedSpec):
    def __init__(self, fid):

        self._fid = fid

        super().__init__()

        self._initialize()

    def _initialize(self):
        pass

    def group_by(self, ):
        pass
    

def read_gene_expression(filename):
    if filename.endswith(".bed"):
        return ParseGeneExpression(io.FileIO(filename))

    if filename.endswith(".bed.gz"):
        return ParseGeneExpression(gzip.GzipFile(filename))

    raise NotImplementedError


class ParseGeneExpression(ExpressionQTLBedSpec):
    def __init__(self, filename):
        pass



def open_predict(filename, mode):

    if mode == "r":
        raise NotImplementedError
    elif mode == "w":
        return WritePredictionBed(io.FileIO(filename, mode))

    raise ValueError("Uninterpretable mode")


class PredictionBedABC(BedSpec):
    """Prediction bed doc string."""
    _req_header_fields = OrderedDict(chrom = str, 
                                     start = int,
                                     end = int,
                                     gene_id = str)


class WritePredictionBed(PredictionBedABC):
    """Buffered write of data to prediction bed file.

    Args:
        fid: (file object)
    """
    def __init__(self, fid):
        super().__init__()

        self._fid = fid
        self._meta_and_header_written = False

        self.meta["data"]="Predicted gene expression"

    def __enter__(self):
        return self

    def __exit__(self, *args):
        if not self.closed:
            self.close()

    @property
    def closed(self):
        return self._fid.closed

    def close(self):
        if not self.closed:
            self._fid.flush()
            self._fid.close()

    def write_meta_data(self, sample_names):
        # instantiate string
        s = ""

        # construct meta data
        for key, val in self.meta.items():
            s += "".join([self._meta_prefix,
                          key,
                          self._meta_data_key_val_delimiter,
                          val,
                          self._new_line])

        # construct header
        s += self._header_prefix
        s += self._field_delimiter.join([*self._req_header_fields.keys(),
                                         *sample_names])
        s = s.encode(encoding=self._encoding)

        if self._fid.write(s) != len(s):
            raise RuntimeError("Bytes written not equal to bytes of string.")

        self._meta_and_header_written = True

    def write_line_record(self, chrom, start, end, name, data):
        """Write line of predictions.

        Args:
            chrom: (str) chromosome name
            start: (int) genomic coordinate of gene beginning
            end: (int) genomic coordinate of gene end
            name: (str) gene id
            data: ((n sample, ) np.ndarray) of floats representing predicted
                gene expression.

        Returns:
            None
        """
        if not self._meta_and_header_written:
            raise ValueError("Write meta data before real data")

        data_str = self._field_delimiter.join([str(w) for w in data])
        chrom = self._new_line + chrom

        record_str = self._field_delimiter.join([chrom, 
                                                str(start), 
                                                str(end), 
                                                name,
                                                data_str])
        record_str = record_str.encode(encoding=self._encoding)

        if self._fid.write(record_str) != len(record_str):
            raise RuntimeError("Bytes written not equal to bytes of string.")

