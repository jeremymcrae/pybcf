''' helpers for writing small BCF files by hand, byte by byte
'''

import struct

from pysam import BGZFile

# bit pattern BCF uses for a missing float, e.g. QUAL of '.'
MISSING_FLOAT = 0x7F800001

def record(pos0=9, n_info=0, n_fmt=0, n_sample=0, chrom=0, rlen=1, qual=MISSING_FLOAT,
           n_allele=1, body=b'\x07\x17A\x00', info=b'', indiv=b''):
    ''' get the bytes of one BCF record

    Args:
        pos0: 0-based position
        n_info: number of INFO fields in `info`
        n_fmt: number of FORMAT fields in `indiv`
        n_sample: number of samples
        chrom: contig index
        rlen: length of the record on the reference
        qual: QUAL, as raw float bits
        n_allele: number of alleles in `body`, including REF
        body: typed ID, alleles and FILTER. The default is a missing ID, REF of
            'A', no ALT, and no FILTER
        info: bytes of the INFO fields
        indiv: bytes of the per-sample FORMAT fields
    '''
    shared = struct.pack('<iiiIII', chrom, pos0, rlen, qual, (n_allele << 16) | n_info,
                         (n_fmt << 24) | n_sample) + body + info
    return struct.pack('<II', len(shared), len(indiv)) + shared + indiv

def write_bcf(path, header, records=()):
    ''' write a bgzip-compressed BCF from VCF header text and record bytes

    Args:
        path: path to write the BCF to
        header: VCF header text, without the terminating NUL
        records: list of record bytes, e.g. from record()
    '''
    text = header.encode('utf8') + b'\x00'
    with BGZFile(str(path), 'wb') as bgz:
        bgz.write(b'BCF\x02\x02' + struct.pack('<I', len(text)) + text)
        bgz.write(b''.join(records))
