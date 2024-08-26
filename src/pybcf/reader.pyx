
import logging
from pathlib import Path

from libcpp cimport bool
from libc.stdint cimport int8_t, uint8_t, uint32_t, int32_t
from libcpp.string cimport string
from libcpp.vector cimport vector

import numpy as np

cdef extern from 'bcf.h' namespace 'bcf':
    cdef cppclass BCF:
        # declare class constructor and methods
        BCF(string path, string index_path) except +
        BCF() except +
        Variant nextvar() except +
        void set_region(string chrom, int32_t start, int32_t end) except +
        Header header

cdef extern from 'header.h' namespace 'bcf':
    cdef cppclass Header:
        # declare class constructor and methods
        Header() except +
        Header(string text) except +
        vector[string] samples
        vector[string] get_contigs()
        vector[string] get_info()
        vector[string] get_filters()
        vector[string] get_formats()

cdef extern from 'info.h' namespace 'bcf':
    cdef struct InfoType:
        int8_t type
        uint32_t offset
    
    cdef cppclass Info:
        # declare class constructor and methods
        Info() except +
        vector[string] get_keys()
        InfoType get_type(string &) except +
        int32_t get_int(uint32_t offset)
        float get_float(uint32_t offset)
        string get_string(uint32_t offset)
        
        vector[int32_t] get_ints(uint32_t offset)
        vector[float] get_floats(uint32_t offset)

cdef extern from 'sample_data.h' namespace 'bcf':
    cdef struct FormatType:
        uint8_t data_type
        uint8_t type_size
        uint32_t offset
        uint32_t n_vals
        bool is_geno
    
    cdef cppclass SampleData:
        # declare class constructor and methods
        SampleData() except +
        vector[string] get_keys()
        FormatType get_type(string &) except +
        vector[int32_t] get_ints(FormatType &)
        vector[float] get_floats(FormatType &)
        vector[string] get_strings(FormatType &)
        uint32_t n_samples
        bool phase_checked
        vector[uint8_t] phase
        vector[uint8_t] missing

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
    cdef list _samples, _contigs, _info, _filters, _formats
    cdef bool samples_done, contigs_done, info_done, filters_done, formats_done
    cdef set_data(self, Header * header):
        ''' assigns Header to python class '''
        self.thisptr = header
        samples_ready = False
    
    @property
    def samples(self):
        ''' list of sample IDs in the BCF
        '''
        if not self.samples_done:
            self._samples = [x.decode('utf8') for x in self.thisptr.samples]
        self.samples_done = True
        return self._samples
    
    @property
    def contigs(self):
        ''' list of possible contigs in the BCF
        '''
        if not self.contigs_done:
            self._contigs = [x.decode('utf8') for x in self.thisptr.get_contigs()]
        self.contigs_done = True
        return self._contigs
    
    @property
    def filters(self):
        ''' list of possible filters in the BCF
        '''
        if not self.filters_done:
            self._filters = [x.decode('utf8') for x in self.thisptr.get_filters()]
        self.filters_done = True
        return self._filters
    
    @property
    def formats(self):
        ''' list of pssible format fields per variant
        '''
        if not self.formats_done:
            self._formats = [x.decode('utf8') for x in self.thisptr.get_formats()]
        self.formats_done = True
        return self._formats
    
    @property
    def info(self):
        ''' list of possible info fields per variant
        '''
        if not self.info_done:
            self._info = [x.decode('utf8') for x in self.thisptr.get_info()]
        self.info_done = True
        return self._info

cdef class BcfInfo:
    cdef Info * thisptr
    cdef set_data(self, Info * info):
        ''' assigns Info to python class '''
        self.thisptr = info
    
    def __contains__(self, key):
        ''' check if key in Info '''
        try:
            self.thisptr.get_type(key.encode('utf8'))
            return True
        except ValueError:
            return False
    
    def keys(self):
        return (x.decode('utf8') for x in self.thisptr.get_keys())
    
    def values(self):
        return (self[x] for x in self.keys())
    
    def __iter__(self):
        ''' get the info keys for a variant
        '''
        for x in self.keys():
            yield x
    
    def __getitem__(self, _key):
        if not self.__contains__(_key):
            raise KeyError(f'unknown INFO field: {_key}')
        
        cdef string key = _key.encode('utf8')
        cdef InfoType info_type = self.thisptr.get_type(key)
        
        cdef int8_t datatype = info_type.type
        cdef uint32_t offset = info_type.offset
        
        if datatype == -1:
            return True
        if datatype == 0:
            return self.thisptr.get_float(offset)
        elif datatype == 1:
            return self.thisptr.get_int(offset)
        elif datatype == 2:
            return self.thisptr.get_string(offset).decode('utf8')
        elif datatype == 3:
            return tuple(self.thisptr.get_floats(offset))
        elif datatype == 4:
            return tuple(self.thisptr.get_ints(offset))

cdef class BcfSampleData:
    cdef SampleData * thisptr
    cdef set_data(self, SampleData * data):
        ''' assigns SampleData to python class '''
        self.thisptr = data
    
    def __contains__(self, key):
        ''' check if key in SampleData '''
        try:
            self.thisptr.get_type(key.encode('utf8'))
            return True
        except ValueError:
            return False
    
    def keys(self):
        return (x.decode('utf8') for x in self.thisptr.get_keys())
    
    def values(self):
        return (self[x] for x in self.keys())
    
    def __iter__(self):
        ''' get the format keys for a variant
        '''
        for x in self.keys():
            yield x
    
    def __getitem__(self, _key):
        if not self.__contains__(_key):
            raise KeyError(f'unknown FORMAT field: {_key}')
        
        cdef string key = _key.encode('utf8')
        
        cdef FormatType fmt_type = self.thisptr.get_type(key)
        cdef uint32_t n_per_sample = fmt_type.n_vals
        cdef uint32_t n_samples = self.thisptr.n_samples
        cdef uint32_t n_vals = n_per_sample * n_samples
        
        cdef vector[int32_t] v_int
        cdef vector[float] v_float
        if fmt_type.data_type == 1 or fmt_type.data_type == 2 or fmt_type.data_type == 3:
            # integer types
            v_int = self.thisptr.get_ints(fmt_type)
            arr = np.asarray(<int32_t [:n_vals]>&v_int[0])
            if n_per_sample > 1:
                arr = np.reshape(arr, (-1, n_per_sample))
            
            mis_val = -1 if key == b'GT' else -(2 ** 31)
            missing = arr == mis_val
            arr = arr.astype(float)
            arr[missing] = float('nan')
            return arr
        elif fmt_type.data_type == 5:
            v_float = self.thisptr.get_floats(fmt_type)
            arr = np.asarray(<float [:n_vals]>&v_float[0])
            if n_per_sample > 1:
                arr = np.reshape(arr, (-1, n_per_sample))
            return arr.copy()
        elif fmt_type.data_type == 7:
            # string type
            return [tuple(x.decode().split(',')) for x in self.thisptr.get_strings(fmt_type)]
        
        raise ValueError(f'unknown datatype: {fmt_type}')
    
    @property
    def phased(self):
        ''' whether genotypes are phased '''
        if not self.thisptr.phase_checked:
            self['GT']
        
        arr = np.asarray(<uint8_t [:self.thisptr.phase.size()]>&self.thisptr.phase[0])
        return arr.astype(np.bool_)

cdef class BcfVariant:
    cdef Variant thisptr
    cdef BcfInfo _info
    cdef BcfSampleData _samples
    def __cinit__(self, BcfReader bcf):
        self.thisptr = bcf.thisptr.nextvar()
        
        # set the _samples attribute
        self._info = BcfInfo()
        self._info.set_data(&self.thisptr.info)
        
        # set the _samples attribute
        self._samples = BcfSampleData()
        self._samples.set_data(&self.thisptr.sample_data)
    
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
        return self.thisptr.varid.decode()
    
    @property
    def qual(self):
        ''' variant quality score (or None)
        '''
        if np.isnan(self.thisptr.qual):
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
        return self._info
    
    @property
    def samples(self):
        ''' variant sample data
        '''
        return self._samples

cdef class BcfReader:
    ''' class to open bcf files from disk, and access variant data within
    '''
    cdef BCF * thisptr
    cdef string path
    cdef bool is_open
    cdef BcfHeader _header
    def __cinit__(self, path, index_path=None):
        path = str(path) if isinstance(path, Path) else path
        index_path = str(index_path) if isinstance(index_path, Path) else index_path
        self.path = path.encode('utf8')
        if index_path is None:
            index_path = b''
        else:
            index_path = index_path.encode('utf8')
        
        logging.debug(f'opening BcfReader from {self.path.decode("utf")}')
        self.thisptr = new BCF(self.path, index_path)
        self.is_open = True
        self._header = BcfHeader()
        self._header.set_data(&self.thisptr.header)
    
    def __close(self):
        if self.is_open:
          del self.thisptr
        self.is_open = False
    
    def __repr__(self):
        return f'BcfReader("{self.path.decode("utf8")}")'
    
    def __iter__(self):
        return self
    
    def __next__(self):
        ''' iterate through all variants in the bcf file
        '''
        try:
            return BcfVariant(self)
        except IndexError:
            raise StopIteration
    
    @property
    def header(self):
      ''' header info from bcf file
      '''
      return self._header
    
    @property
    def samples(self):
      ''' list of samples in the bcf file
      '''
      return self.header.samples
    
    def fetch(self, chrom, int32_t start=1, int32_t stop=2**29):
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
        
        try:
            self.thisptr.set_region(chrom.encode('utf8'), start, stop)
        except IndexError:
            pass
        
        return iter(self)
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.__close()
        return False
