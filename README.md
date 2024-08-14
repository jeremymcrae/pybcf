
### pybcf

This is a package for reading BCF files (binary call format). Similar to pysam/cyvcf2, 
but limited to reading from BCF files only, and oriented around getting genotype 
data into numpy arrays quickly. 

This performs relatively better at larger sample sizes e.g. it's 7X faster than
pysam with 2k samples, but improves to 50X faster than pysam with 30k samples. 
Currently it's slower than pysam for BCFs without sample level data, but BCFs
without genotypes aren't generally a limiting factor e.g. it still parses 
300k variants/second under those conditions.

```py
from pybcf import BcfReader

bcf = BcfReader(bcf_path)

sample_ids = bcf.samples
contigs = bcf.header.contigs

for var in bcf:
    # the usual attributes are available e.g.
    # var.chrom, var.pos, var.ref, var.alts, var.info['AF']
    
    # sample data is accessed as numpy arrays via the format keys
    keys = list(var.samples)
    genotypes = var.samples['GT']  # as n x 2 numpy array, missing=nan

# or fetch from random regions if bcf is indexed
for var in bcf.fetch('chr1', 10000, 200000):
    print(var.chrom, var.pos)
```

### Limitations
 - doesn't work with uncompressed BCFs
 - extracting info fields is a little slow
