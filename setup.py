
import glob
import io
from setuptools import setup
import sys
import os

from distutils.core import Extension
from distutils.ccompiler import new_compiler
from Cython.Build import cythonize

EXTRA_COMPILE_ARGS = ['-std=c++11', '-I/usr/include']
EXTRA_LINK_ARGS = []
if sys.platform == "darwin":
    EXTRA_COMPILE_ARGS += [
        "-stdlib=libc++",
        "-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1",
        "-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include",
        ]
    EXTRA_LINK_ARGS += [
        "-L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib",
        ]

def flatten(*lists):
    return [str(x) for sublist in lists for x in sublist]

def build_zlib():
    ''' compile zlib code to object files for linking with bgen on windows
    
    Returns:
        list of paths to compiled object code
    '''
    include_dirs = ['src/zlib/']
    sources = list(glob.glob('src/zlib/*.c'))
    extra_compile_args = ['/O2']
    
    cc = new_compiler()
    return cc.compile(sources, include_dirs=include_dirs,
        extra_preargs=extra_compile_args)

def get_gzstream_path():
    ''' workaround for building gzstream on windows

    cython on windows didn't like the .C extension for gzstream. This just
    renames the file (on windows only), and returns the relative path.
    '''
    gzstream_path = 'src/gzstream/gzstream.C'
    if sys.platform == 'win32':
        gzstream_win_path = 'src/gzstream/gzstream.cpp'
        try:
            os.rename(gzstream_path, gzstream_win_path)
        except FileNotFoundError:
            pass  # avoid error on github actions
        gzstream_path = gzstream_win_path
    return gzstream_path

def scrub_gzstream():
    ''' workaround for compilation error on macos
    
    compiling gzstream requires the corresponding gzstream.h file, but if we 
    include the gzstream directory in the include dirs, then clang complains
    about the version file in the gzstream folder. If we remove the gzstream
    directory from the include dirs, then clang complains about the missing
    gzstream.h. This is because gzstream.C identifies it's header file with
    angle brackets. Replacing the angle brackets in that line seems to work.
    '''
    with open(get_gzstream_path(), 'rt') as handle:
        lines = handle.readlines()
    
    with open(get_gzstream_path(), 'wt') as handle:
        for line in lines:
            if line == '#include <gzstream.h>\n':
                line = '#include "gzstream.h"\n'
            handle.write(line)

if sys.platform == 'win32':
    zlib, libs = build_zlib(), []
else:
    zlib, libs = [], ['z']

scrub_gzstream()

ext = cythonize([
    Extension('pybcf.reader',
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        sources=['src/pybcf/reader.pyx',
            get_gzstream_path(),
            'src/bcf.cpp',
            'src/header.cpp',
            'src/info.cpp',
            'src/sample_data.cpp',
            'src/variant.cpp'],
        extra_objects=zlib,
        include_dirs=['src', 'src/zlib'],
        libraries=libs,
        language='c++'),
    ])

setup(name='pybcf',
    description='Package for loading data from bcf files',
    long_description=io.open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    version='1.0.0',
    author='Jeremy McRae',
    author_email='jmcrae@illumina.com',
    license="MIT",
    url='https://github.com/jeremymcrae/pybcf',
    packages=['pybcf'],
    package_dir={'': 'src'},
    install_requires=[
        'numpy',
    ],
    extras_require={
        'test': [
            'pysam',
         ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    ext_modules=ext,
    test_loader='unittest:TestLoader',
    )