import logging
from pathlib import Path

from libcpp cimport bool
from libc.math cimport isnan
from libc.stddef cimport size_t
from libc.stdint cimport uint8_t, uint32_t, int32_t, uint64_t, INT32_MIN, INT32_MAX
from libcpp.string cimport string
from libcpp.vector cimport vector

import numpy as np

cdef extern from 'types.h' namespace 'bcf':
    cdef enum Types:
        flag
        int8
        int16
        int32
        float_
        char_

cdef extern from 'bcf.h' namespace 'bcf':
    cdef cppclass BCF:
        # declare class constructor and methods
        BCF(string path, string index_path) except +
        BCF() except +
        bool nextvar(Variant &) except +
        void set_region(string chrom, int32_t start, int32_t end) except +
        void close() except +
        Header header

cdef extern from 'header.h' namespace 'bcf':
    cdef cppclass Header:
        # declare class constructor and methods
        Header() except +
        Header(const string & text) except +
        vector[string] samples
        vector[string] get_contigs() except +
        vector[string] get_info() except +
        vector[string] get_filters() except +
        vector[string] get_formats() except +

cdef extern from 'info.h' namespace 'bcf':
    cdef enum InfoKind:
        info_missing
        info_flag
        info_scalar_float
        info_scalar_int
        info_string
        info_vector_float
        info_vector_int
        info_vector_string

    cdef struct InfoType:
        InfoKind type
        uint8_t type_size
        uint32_t offset
        uint32_t n_vals

    cdef cppclass Info:
        # declare class constructor and methods
        Info() except +
        size_t size() except +
        vector[string] get_keys() except +
        bool get_type(const string &, InfoType &) except +
        int32_t get_int(const InfoType &) except +
        float get_float(const InfoType &) except +
        string get_string(const InfoType &) except +

        vector[int32_t] get_ints(const InfoType &) except +
        vector[float] get_floats(const InfoType &) except +

cdef extern from 'sample_data.h' namespace 'bcf':
    cdef struct FormatType:
        Types data_type
        uint8_t type_size
        uint32_t offset
        uint32_t n_vals
        bool is_geno

    cdef cppclass SampleData:
        # declare class constructor and methods
        SampleData() except +
        size_t size() except +
        vector[string] get_keys() except +
        bool get_type(const string &, FormatType &) except +
        vector[int32_t] get_ints(const FormatType &) except +
        vector[float] get_floats(const FormatType &) except +
        vector[string] get_strings(const FormatType &) except +
        uint32_t n_samples
        bool phase_checked
        vector[uint8_t] phase

cdef extern from 'variant.h' namespace 'bcf':
    cdef cppclass Variant:
        # declare class constructor and methods
        Variant() except +
        string chrom
        int pos
        string ref
        vector[string] alts
        float qual
        string varid
        vector[string] filters
        Info info
        SampleData sample_data

cdef class BcfHeader:
    cdef Header * thisptr
    cdef BcfReader _reader
    cdef list _samples, _contigs, _info, _filters, _formats
    cdef set_data(self, Header * header, BcfReader reader):
        ''' assigns Header to python class '''
        self.thisptr = header
        self._reader = reader

    @property
    def samples(self):
        ''' list of sample IDs in the BCF
        '''
        if self._samples is None:
            self._samples = [x.decode('utf8') for x in self.thisptr.samples]
        return self._samples

    @property
    def contigs(self):
        ''' list of possible contigs in the BCF
        '''
        if self._contigs is None:
            self._contigs = [x.decode('utf8') for x in self.thisptr.get_contigs()]
        return self._contigs

    @property
    def filters(self):
        ''' list of possible filters in the BCF
        '''
        if self._filters is None:
            self._filters = [x.decode('utf8') for x in self.thisptr.get_filters()]
        return self._filters

    @property
    def formats(self):
        ''' list of possible format fields per variant
        '''
        if self._formats is None:
            self._formats = [x.decode('utf8') for x in self.thisptr.get_formats()]
        return self._formats

    @property
    def info(self):
        ''' list of possible info fields per variant
        '''
        if self._info is None:
            self._info = [x.decode('utf8') for x in self.thisptr.get_info()]
        return self._info

cdef class BcfInfo:
    cdef Info * thisptr
    cdef BcfVariant _variant
    cdef set_data(self, Info * info, BcfVariant variant):
        ''' assigns Info to python class '''
        self.thisptr = info
        self._variant = variant

    def __len__(self):
        return self.thisptr.size()

    def __contains__(self, key):
        ''' check if key in Info '''
        cdef InfoType info_type
        return self.thisptr.get_type(key.encode('utf8'), info_type)

    def keys(self):
        return (x.decode('utf8') for x in self.thisptr.get_keys())

    def values(self):
        return (self[x] for x in self.keys())

    def items(self):
        return ((x, self[x]) for x in self.keys())

    def __iter__(self):
        ''' get the info keys for a variant
        '''
        return self.keys()

    cdef _decode_value(self, InfoType info_type):
        cdef InfoKind datatype = info_type.type
        cdef float f_val
        cdef int32_t i_val

        if datatype == info_missing:
            return None
        elif datatype == info_flag:
            return True
        elif datatype == info_scalar_float:
            f_val = self.thisptr.get_float(info_type)
            return None if isnan(f_val) else f_val
        elif datatype == info_scalar_int:
            i_val = self.thisptr.get_int(info_type)
            return None if i_val == INT32_MIN else i_val
        elif datatype == info_string:
            return self.thisptr.get_string(info_type).decode('utf8')
        elif datatype == info_vector_string:
            return tuple(self.thisptr.get_string(info_type).decode('utf8').split(','))
        elif datatype == info_vector_float:
            return tuple(None if isnan(x) else x for x in self.thisptr.get_floats(info_type))
        elif datatype == info_vector_int:
            return tuple(None if x == INT32_MIN else x for x in self.thisptr.get_ints(info_type))

    def get(self, key, default=None):
        cdef InfoType info_type
        if not self.thisptr.get_type(key.encode('utf8'), info_type):
            return default
        return self._decode_value(info_type)

    def __getitem__(self, _key):
        cdef InfoType info_type
        if not self.thisptr.get_type(_key.encode('utf8'), info_type):
            raise KeyError(f'unknown INFO field: {_key}')
        return self._decode_value(info_type)

cdef class BcfSampleData:
    cdef SampleData * thisptr
    cdef BcfVariant _variant
    cdef set_data(self, SampleData * data, BcfVariant variant):
        ''' assigns SampleData to python class '''
        self.thisptr = data
        self._variant = variant

    def __len__(self):
        return self.thisptr.size()

    def __contains__(self, key):
        ''' check if key in SampleData '''
        cdef FormatType fmt_type
        return self.thisptr.get_type(key.encode('utf8'), fmt_type)

    def keys(self):
        return (x.decode('utf8') for x in self.thisptr.get_keys())

    def values(self):
        return (self[x] for x in self.keys())

    def items(self):
        return ((x, self[x]) for x in self.keys())

    def __iter__(self):
        ''' get the format keys for a variant
        '''
        return self.keys()

    cdef _decode_value(self, FormatType fmt_type):
        cdef uint32_t n_per_sample = fmt_type.n_vals
        cdef uint32_t n_samples = self.thisptr.n_samples
        cdef uint64_t n_vals = <uint64_t>n_per_sample * <uint64_t>n_samples
        cdef tuple shape
        if n_per_sample > 1:
            shape = (n_samples, n_per_sample)
        else:
            shape = (n_samples,)

        cdef vector[int32_t] v_int
        cdef vector[float] v_float
        if fmt_type.data_type in (int8, int16, int32):
            if n_vals == 0:
                return np.full(shape, np.nan, dtype=float)
            v_int = self.thisptr.get_ints(fmt_type)
            arr = np.asarray(<int32_t [:n_vals]>&v_int[0]).reshape(shape)

            mis_val = -1 if fmt_type.is_geno else INT32_MIN
            missing = arr == mis_val
            arr = arr.astype(float)
            arr[missing] = float('nan')
            return arr
        elif fmt_type.data_type == float_:
            if n_vals == 0:
                return np.full(shape, np.nan, dtype=np.float32)
            v_float = self.thisptr.get_floats(fmt_type)
            return np.asarray(<float [:n_vals]>&v_float[0]).reshape(shape).copy()
        elif fmt_type.data_type == char_:
            return [tuple(x.decode('utf8').split(',')) for x in self.thisptr.get_strings(fmt_type)]

        raise ValueError(f'unknown datatype: {<int>fmt_type.data_type}')

    def get(self, key, default=None):
        cdef FormatType fmt_type
        if not self.thisptr.get_type(key.encode('utf8'), fmt_type):
            return default
        return self._decode_value(fmt_type)

    def __getitem__(self, _key):
        cdef FormatType fmt_type
        if not self.thisptr.get_type(_key.encode('utf8'), fmt_type):
            raise KeyError(f'unknown FORMAT field: {_key}')
        return self._decode_value(fmt_type)

    @property
    def phased(self):
        ''' whether genotypes are phased '''
        if self.thisptr.n_samples == 0:
            return np.empty(0, dtype=np.bool_)
        if not self.thisptr.phase_checked:
            self.get('GT')
        if self.thisptr.phase.empty():
            return np.empty(0, dtype=np.bool_)

        arr = np.asarray(<uint8_t [:self.thisptr.phase.size()]>&self.thisptr.phase[0])
        return arr.astype(np.bool_)

cdef class BcfVariant:
    cdef Variant thisptr
    cdef BcfReader _reader
    cdef BcfInfo _info
    cdef BcfSampleData _samples
    def __cinit__(self, BcfReader bcf):
        self._reader = bcf

    @property
    def chrom(self):
        ''' chromosome/contig name
        '''
        return self.thisptr.chrom.decode('utf8')

    @property
    def contig(self):
        ''' alias for .chrom attribute
        '''
        return self.chrom

    @property
    def pos(self):
        ''' variant nucleotide position
        '''
        return self.thisptr.pos

    @property
    def ref(self):
        ''' reference allele
        '''
        return self.thisptr.ref.decode('utf8')

    @property
    def alts(self):
        ''' list of alternate alleles
        '''
        if self.thisptr.alts.size() == 0:
            return None
        return tuple(x.decode('utf8') for x in self.thisptr.alts)

    @property
    def id(self):
        ''' variant ID
        '''
        if self.thisptr.varid.size() == 0:
            return None
        return self.thisptr.varid.decode('utf8')

    @property
    def qual(self):
        ''' variant quality score (or None)
        '''
        if isnan(self.thisptr.qual):
            return None
        return self.thisptr.qual

    @property
    def filter(self):
        ''' list of filter information
        '''
        return [x.decode('utf8') for x in self.thisptr.filters]

    @property
    def info(self):
        ''' variant info data
        '''
        if self._info is None:
            self._info = BcfInfo()
            self._info.set_data(&self.thisptr.info, self)
        return self._info

    @property
    def samples(self):
        ''' variant sample data
        '''
        if self._samples is None:
            self._samples = BcfSampleData()
            self._samples.set_data(&self.thisptr.sample_data, self)
        return self._samples

cdef class BcfReader:
    ''' class to open bcf files from disk, and access variant data within
    '''
    cdef BCF * thisptr
    cdef string path
    cdef bool is_open
    cdef BcfHeader _header
    def __cinit__(self, path, index_path=None):
        self.thisptr = NULL
        self.is_open = False
        path = str(path) if isinstance(path, Path) else path
        index_path = str(index_path) if isinstance(index_path, Path) else index_path
        self.path = path.encode('utf8')
        if index_path is None:
            index_path = b''
        else:
            index_path = index_path.encode('utf8')

        logging.debug(f'opening BcfReader from {self.path.decode("utf8")}')
        self.thisptr = new BCF(self.path, index_path)
        self.is_open = True

    def __dealloc__(self):
        if self.thisptr != NULL:
            del self.thisptr
            self.thisptr = NULL
        self.is_open = False

    def close(self):
        ''' close the bcf file '''
        if self.is_open and self.thisptr != NULL:
            self.thisptr.close()
        self.is_open = False

    def __repr__(self):
        return f'BcfReader("{self.path.decode("utf8")}")'

    def __iter__(self):
        return self

    def __next__(self):
        ''' iterate through all variants in the bcf file
        '''
        if not self.is_open:
            raise ValueError('bcf is closed')
        cdef BcfVariant var = BcfVariant(self)
        if not self.thisptr.nextvar(var.thisptr):
            raise StopIteration
        return var

    @property
    def header(self):
        ''' header info from bcf file
        '''
        if self._header is None:
            self._header = BcfHeader()
            self._header.set_data(&self.thisptr.header, self)
        return self._header

    @property
    def samples(self):
        ''' list of samples in the bcf file
        '''
        return self.header.samples

    def fetch(self, chrom, start=None, stop=None):
        ''' fetches all variants within a genomic region

        Args:
            chrom: chromosome that variants must be on
            start: start nucleotide of region. If None, starts with first variant on chrom
            stop: end nucleotide of region. If None, continues to last variant on chrom

        Yields:
            BcfVars for variants within the genome region
        '''
        if not self.is_open:
            raise ValueError('bcf is closed')

        cdef int32_t region_start = 1 if start is None else start
        cdef int32_t region_stop = INT32_MAX if stop is None else stop
        self.thisptr.set_region(chrom.encode('utf8'), region_start, region_stop)
        return iter(self)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False
