# cython: boundscheck=False, wraparound=False
import logging
from pathlib import Path

from cpython.unicode cimport PyUnicode_AsUTF8AndSize
from libcpp cimport bool
from libc.math cimport isnan
from libc.stddef cimport size_t
from libc.stdint cimport uint8_t, uint32_t, int32_t, uint64_t, INT32_MIN, INT32_MAX
from libc.string cimport memcpy
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
        BCF(const string & path, const string & index_path) except +
        BCF() except +
        bool nextvar(Variant &) except +
        void set_region(const string & chrom, int32_t start, int32_t end) except +
        void close() except +
        Header header

cdef extern from 'header.h' namespace 'bcf':
    cdef cppclass Header:
        Header() except +
        Header(const string & text) except +
        vector[string] contigs
        vector[string] filters
        vector[string] samples
        vector[string] get_contigs() except +
        vector[string] get_info() except +
        vector[string] get_filters() except +
        vector[string] get_formats() except +
        vector[string] get_info_by_idx() except +
        vector[string] get_formats_by_idx() except +

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
        Info() except +
        size_t size() except +
        uint32_t key_id_at(size_t) except +
        const InfoType & type_at(size_t) except +
        bool get_type(const char *, size_t, InfoType &) except +
        int32_t get_int(const InfoType &) except +
        int32_t get_raw_int(const InfoType &) except +
        float get_float(const InfoType &) except +
        uint32_t get_raw_float_bits(const InfoType &) except +
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
        SampleData() except +
        size_t size() except +
        uint32_t key_id_at(size_t) except +
        const FormatType & type_at(size_t) except +
        bool get_type(const char *, size_t, FormatType &) except +
        void read_ints(const FormatType &, double *) except +
        void read_floats(const FormatType &, float *) except +
        void read_phase(const FormatType &, uint8_t *) except +
        vector[string] get_strings(const FormatType &) except +
        uint32_t n_samples

cdef extern from 'variant.h' namespace 'bcf':
    cdef cppclass Variant:
        Variant() except +
        uint32_t chrom_id
        int32_t pos
        string ref
        vector[string] alts
        float qual
        string varid
        vector[uint32_t] filter_ids
        Info info
        SampleData sample_data

cdef inline tuple _decode_string_vec(const vector[string] & vec):
    cdef size_t i, n = vec.size()
    return tuple([vec[i].decode('utf8') if not vec[i].empty() else None for i in range(n)])

cdef class BcfHeader:
    cdef Header * thisptr
    cdef BcfReader _reader
    cdef list _samples, _contigs, _info, _filters, _formats

    def __init__(self, *args, **kwargs):
        raise TypeError('cannot instantiate BcfHeader directly')

    cdef set_data(self, Header * header, BcfReader reader):
        ''' assigns Header to python class '''
        self.thisptr = header
        self._reader = reader

    @property
    def samples(self):
        ''' list of sample IDs in the BCF
        '''
        cdef size_t i, n
        if self._samples is None:
            n = self.thisptr.samples.size()
            self._samples = [self.thisptr.samples[i].decode('utf8') for i in range(n)]
        return self._samples

    @property
    def contigs(self):
        ''' list of possible contigs in the BCF
        '''
        cdef vector[string] vals
        cdef size_t i, n
        if self._contigs is None:
            vals = self.thisptr.get_contigs()
            n = vals.size()
            self._contigs = [vals[i].decode('utf8') for i in range(n)]
        return self._contigs

    @property
    def filters(self):
        ''' list of possible filters in the BCF
        '''
        cdef vector[string] vals
        cdef size_t i, n
        if self._filters is None:
            vals = self.thisptr.get_filters()
            n = vals.size()
            self._filters = [vals[i].decode('utf8') for i in range(n)]
        return self._filters

    @property
    def formats(self):
        ''' list of possible format fields per variant
        '''
        cdef vector[string] vals
        cdef size_t i, n
        if self._formats is None:
            vals = self.thisptr.get_formats()
            n = vals.size()
            self._formats = [vals[i].decode('utf8') for i in range(n)]
        return self._formats

    @property
    def info(self):
        ''' list of possible info fields per variant
        '''
        cdef vector[string] vals
        cdef size_t i, n
        if self._info is None:
            vals = self.thisptr.get_info()
            n = vals.size()
            self._info = [vals[i].decode('utf8') for i in range(n)]
        return self._info

cdef class BcfInfo:
    cdef Info * thisptr
    cdef BcfVariant _variant

    def __init__(self, *args, **kwargs):
        raise TypeError('cannot instantiate BcfInfo directly')

    cdef set_data(self, Info * info, BcfVariant variant):
        ''' assigns Info to python class '''
        self.thisptr = info
        self._variant = variant

    def __len__(self):
        return self.thisptr.size()

    def __contains__(self, key):
        ''' check if key in Info '''
        if not isinstance(key, str):
            return False
        cdef Py_ssize_t key_len = 0
        cdef const char * key_ptr = PyUnicode_AsUTF8AndSize(key, &key_len)
        cdef InfoType info_type
        return self.thisptr.get_type(key_ptr, <size_t>key_len, info_type)

    def keys(self):
        cdef size_t i, n = self.thisptr.size()
        cdef tuple names = self._variant._reader._get_info_lookup()
        return [names[self.thisptr.key_id_at(i)] for i in range(n)]

    def values(self):
        cdef size_t i, n = self.thisptr.size()
        return [self._decode_value(self.thisptr.type_at(i)) for i in range(n)]

    def items(self):
        cdef size_t i, n = self.thisptr.size()
        cdef tuple names = self._variant._reader._get_info_lookup()
        return [(names[self.thisptr.key_id_at(i)],
                 self._decode_value(self.thisptr.type_at(i))) for i in range(n)]

    def __iter__(self):
        ''' get the info keys for a variant
        '''
        return iter(self.keys())

    cdef _decode_value(self, InfoType info_type):
        cdef InfoKind datatype = info_type.type
        cdef float f_val
        cdef uint32_t f_bits
        cdef int32_t i_val
        cdef vector[float] v_float
        cdef vector[int32_t] v_int
        cdef size_t i, n

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
            if info_type.n_vals == 0:
                return ()
            if info_type.n_vals == 1:
                f_bits = self.thisptr.get_raw_float_bits(info_type)
                if f_bits == 0x7f800002:
                    return ()
                if f_bits == 0x7f800001:
                    return (None,)
                memcpy(&f_val, &f_bits, 4)
                return (f_val,)
            v_float = self.thisptr.get_floats(info_type)
            n = v_float.size()
            return tuple([None if isnan(v_float[i]) else v_float[i] for i in range(n)])
        elif datatype == info_vector_int:
            if info_type.n_vals == 0:
                return ()
            if info_type.n_vals == 1:
                i_val = self.thisptr.get_raw_int(info_type)
                if i_val == INT32_MIN + 1:
                    return ()
                return (None if i_val == INT32_MIN else i_val,)
            v_int = self.thisptr.get_ints(info_type)
            n = v_int.size()
            return tuple([None if v_int[i] == INT32_MIN else v_int[i] for i in range(n)])

    def get(self, key, default=None):
        if not isinstance(key, str):
            return default
        cdef Py_ssize_t key_len = 0
        cdef const char * key_ptr = PyUnicode_AsUTF8AndSize(key, &key_len)
        cdef InfoType info_type
        if not self.thisptr.get_type(key_ptr, <size_t>key_len, info_type):
            return default
        return self._decode_value(info_type)

    def __getitem__(self, _key):
        if not isinstance(_key, str):
            raise KeyError(f'unknown INFO field: {_key}')
        cdef Py_ssize_t key_len = 0
        cdef const char * key_ptr = PyUnicode_AsUTF8AndSize(_key, &key_len)
        cdef InfoType info_type
        if not self.thisptr.get_type(key_ptr, <size_t>key_len, info_type):
            raise KeyError(f'unknown INFO field: {_key}')
        return self._decode_value(info_type)

cdef class BcfSampleData:
    cdef SampleData * thisptr
    cdef BcfVariant _variant

    def __init__(self, *args, **kwargs):
        raise TypeError('cannot instantiate BcfSampleData directly')

    cdef set_data(self, SampleData * data, BcfVariant variant):
        ''' assigns SampleData to python class '''
        self.thisptr = data
        self._variant = variant

    def __len__(self):
        return self.thisptr.size()

    def __contains__(self, key):
        ''' check if key in SampleData '''
        if not isinstance(key, str):
            return False
        cdef Py_ssize_t key_len = 0
        cdef const char * key_ptr = PyUnicode_AsUTF8AndSize(key, &key_len)
        cdef FormatType fmt_type
        return self.thisptr.get_type(key_ptr, <size_t>key_len, fmt_type)

    def keys(self):
        cdef size_t i, n = self.thisptr.size()
        cdef tuple names = self._variant._reader._get_format_lookup()
        return [names[self.thisptr.key_id_at(i)] for i in range(n)]

    def values(self):
        cdef size_t i, n = self.thisptr.size()
        return [self._decode_value(self.thisptr.type_at(i)) for i in range(n)]

    def items(self):
        cdef size_t i, n = self.thisptr.size()
        cdef tuple names = self._variant._reader._get_format_lookup()
        return [(names[self.thisptr.key_id_at(i)],
                 self._decode_value(self.thisptr.type_at(i))) for i in range(n)]

    def __iter__(self):
        ''' get the format keys for a variant
        '''
        return iter(self.keys())

    cdef _decode_value(self, FormatType fmt_type):
        cdef uint32_t n_per_sample = fmt_type.n_vals
        cdef uint32_t n_samples = self.thisptr.n_samples
        cdef uint64_t n_vals = <uint64_t>n_per_sample * <uint64_t>n_samples
        cdef tuple shape
        cdef double[::1] i_arr
        cdef float[::1] f_arr
        cdef vector[string] v_str
        cdef size_t i, n

        if fmt_type.data_type == int8 or fmt_type.data_type == int16 or fmt_type.data_type == int32:
            if n_vals == 0:
                if n_per_sample > 1:
                    shape = (n_samples, n_per_sample)
                else:
                    shape = (n_samples,)
                return np.full(shape, np.nan, dtype=np.float64)
            raw_arr = np.empty(n_vals, dtype=np.float64, order='C')
            i_arr = raw_arr
            self.thisptr.read_ints(fmt_type, &i_arr[0])
            return raw_arr.reshape((n_samples, n_per_sample)) if n_per_sample > 1 else raw_arr
        elif fmt_type.data_type == float_:
            if n_vals == 0:
                if n_per_sample > 1:
                    shape = (n_samples, n_per_sample)
                else:
                    shape = (n_samples,)
                return np.full(shape, np.nan, dtype=np.float32)
            raw_arr = np.empty(n_vals, dtype=np.float32, order='C')
            f_arr = raw_arr
            self.thisptr.read_floats(fmt_type, &f_arr[0])
            return raw_arr.reshape((n_samples, n_per_sample)) if n_per_sample > 1 else raw_arr
        elif fmt_type.data_type == char_:
            v_str = self.thisptr.get_strings(fmt_type)
            n = v_str.size()
            return [tuple(v_str[i].decode('utf8').split(',')) for i in range(n)]

        raise ValueError(f'unknown datatype: {<int>fmt_type.data_type}')

    def get(self, key, default=None):
        if not isinstance(key, str):
            return default
        cdef Py_ssize_t key_len = 0
        cdef const char * key_ptr = PyUnicode_AsUTF8AndSize(key, &key_len)
        cdef FormatType fmt_type
        if not self.thisptr.get_type(key_ptr, <size_t>key_len, fmt_type):
            return default
        return self._decode_value(fmt_type)

    def __getitem__(self, _key):
        if not isinstance(_key, str):
            raise KeyError(f'unknown FORMAT field: {_key}')
        cdef Py_ssize_t key_len = 0
        cdef const char * key_ptr = PyUnicode_AsUTF8AndSize(_key, &key_len)
        cdef FormatType fmt_type
        if not self.thisptr.get_type(key_ptr, <size_t>key_len, fmt_type):
            raise KeyError(f'unknown FORMAT field: {_key}')
        return self._decode_value(fmt_type)

    @property
    def phased(self):
        ''' whether genotypes are phased '''
        if self.thisptr.n_samples == 0:
            return np.empty(0, dtype=np.bool_)
        cdef FormatType fmt_type
        if not self.thisptr.get_type(b'GT', 2, fmt_type):
            return np.empty(0, dtype=np.bool_)
        raw_arr = np.empty(self.thisptr.n_samples, dtype=np.uint8, order='C')
        cdef uint8_t[::1] p_arr = raw_arr
        self.thisptr.read_phase(fmt_type, &p_arr[0])
        return raw_arr.view(np.bool_)

cdef class BcfVariant:
    cdef Variant thisptr
    cdef BcfReader _reader
    cdef BcfInfo _info
    cdef BcfSampleData _samples
    def __cinit__(self, BcfReader bcf):
        self._reader = bcf

    def __init__(self, *args, **kwargs):
        raise TypeError('cannot instantiate BcfVariant directly')

    @property
    def chrom(self):
        ''' chromosome/contig name
        '''
        return self._reader._get_contig_lookup()[self.thisptr.chrom_id]

    @property
    def contig(self):
        ''' alias for .chrom attribute
        '''
        return self._reader._get_contig_lookup()[self.thisptr.chrom_id]

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
        cdef size_t i, n = self.thisptr.alts.size()
        if n == 0:
            return None
        if n == 1:
            return (self.thisptr.alts[0].decode('utf8'),)
        return tuple([self.thisptr.alts[i].decode('utf8') for i in range(n)])

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
        cdef size_t i, n = self.thisptr.filter_ids.size()
        if n == 0:
            return []
        cdef tuple names = self._reader._get_filter_lookup()
        if n == 1:
            return [names[self.thisptr.filter_ids[0]]]
        return [names[self.thisptr.filter_ids[i]] for i in range(n)]

    @property
    def info(self):
        ''' variant info data
        '''
        if self._info is None:
            self._info = BcfInfo.__new__(BcfInfo)
            self._info.set_data(&self.thisptr.info, self)
        return self._info

    @property
    def samples(self):
        ''' variant sample data
        '''
        if self._samples is None:
            self._samples = BcfSampleData.__new__(BcfSampleData)
            self._samples.set_data(&self.thisptr.sample_data, self)
        return self._samples

cdef class BcfReader:
    ''' class to open bcf files from disk, and access variant data within
    '''
    cdef BCF * thisptr
    cdef string path
    cdef bool is_open
    cdef BcfHeader _header
    cdef tuple _contig_lookup
    cdef tuple _filter_lookup
    cdef tuple _info_lookup
    cdef tuple _format_lookup
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

    cdef inline tuple _get_contig_lookup(self):
        if self._contig_lookup is None:
            self._contig_lookup = _decode_string_vec(self.thisptr.header.contigs)
        return self._contig_lookup

    cdef inline tuple _get_filter_lookup(self):
        if self._filter_lookup is None:
            self._filter_lookup = _decode_string_vec(self.thisptr.header.filters)
        return self._filter_lookup

    cdef inline tuple _get_info_lookup(self):
        if self._info_lookup is None:
            self._info_lookup = _decode_string_vec(self.thisptr.header.get_info_by_idx())
        return self._info_lookup

    cdef inline tuple _get_format_lookup(self):
        if self._format_lookup is None:
            self._format_lookup = _decode_string_vec(self.thisptr.header.get_formats_by_idx())
        return self._format_lookup

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
        cdef BcfVariant var = BcfVariant.__new__(BcfVariant, self)
        if not self.thisptr.nextvar(var.thisptr):
            raise StopIteration
        return var

    @property
    def header(self):
        ''' header info from bcf file
        '''
        if self._header is None:
            self._header = BcfHeader.__new__(BcfHeader)
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
