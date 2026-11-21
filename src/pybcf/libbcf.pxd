# Cython declarations for the C++ BCF reader classes in `src/*.h`. These only tell Cython
# which parts of the C++ classes `reader.pyx` uses; the classes themselves are in C++.
from libcpp cimport bool
from libc.stddef cimport size_t
from libc.stdint cimport uint8_t, uint32_t, int32_t
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from 'types.h' namespace 'bcf':
    cdef uint32_t BCF_FLOAT_MISSING
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
        bool nextvar(Variant &) except +
        void set_region(const string & chrom, int32_t start, int32_t end) except +
        void close() except +
        Header header

cdef extern from 'header.h' namespace 'bcf':
    cdef struct InfoField:
        string id

    cdef struct FormatField:
        string id

    cdef cppclass Header:
        vector[string] contigs
        vector[InfoField] info
        vector[string] filters
        vector[FormatField] format
        vector[string] samples

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

    cdef struct InfoLayout:
        InfoKind kind
        uint8_t type_size
        uint32_t offset
        uint32_t n_vals

    cdef struct NamedInfo:
        uint32_t id_idx
        InfoLayout layout

    cdef cppclass Info:
        size_t size() except +
        vector[NamedInfo] & get_fields() except +
        bool find(const char *, size_t, InfoLayout &) except +
        int32_t get_int(const InfoLayout &) except +
        uint32_t get_float_bits(const InfoLayout &) except +
        string get_string(const InfoLayout &) except +
        vector[int32_t] get_ints(const InfoLayout &) except +
        vector[uint32_t] get_floats_bits(const InfoLayout &) except +

cdef extern from 'sample_data.h' namespace 'bcf':
    cdef struct FormatLayout:
        Types data_type
        uint8_t type_size
        uint32_t offset
        uint32_t n_vals
        bool is_gt

    cdef struct NamedFormat:
        uint32_t id_idx
        FormatLayout layout

    cdef cppclass SampleData:
        size_t size() except +
        vector[NamedFormat] & get_fields() except +
        bool find(const char *, size_t, FormatLayout &) except +
        void read_ints(const FormatLayout &, double *) except +
        void read_floats(const FormatLayout &, float *) except +
        void read_phase(const FormatLayout &, uint8_t *) except +
        vector[string] get_strings(const FormatLayout &) except +
        uint32_t n_samples

cdef extern from 'variant.h' namespace 'bcf':
    cdef cppclass Variant:
        Variant() except +
        uint32_t chrom_id
        int32_t pos
        string ref
        vector[string] alts
        uint32_t qual_bits
        string varid
        vector[uint32_t] filter_ids
        Info info
        SampleData sample_data
