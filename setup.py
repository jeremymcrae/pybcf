
import io
from pathlib import Path
from setuptools import setup
import subprocess
import sys
import os
import platform

from distutils.core import Extension
from Cython.Build import cythonize

EXTRA_LINK_ARGS = []
EXTRA_COMPILE_ARGS = []
if sys.platform == 'linux':
    EXTRA_COMPILE_ARGS = ['-std=c++11', '-I/usr/include']
elif sys.platform == "darwin":
    EXTRA_COMPILE_ARGS += [
        "-stdlib=libc++",
        "-std=c++11",
        '-mavx',   # required for macosx testing on github actions
        '-mavx2',  # required for macosx testing on github actions
        "-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1",
        "-I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include",
        ]
    EXTRA_LINK_ARGS += [
        "-L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib",
        ]
elif sys.platform == "win32":
    EXTRA_COMPILE_ARGS += ['/std:c++14', '/O2']

if platform.machine() == 'x86_64':
    EXTRA_COMPILE_ARGS += ['-mavx', '-mavx2']

def build_zlib():
    ''' compile zlib code to object files for linking
    
    Returns:
        list of paths to compiled object code
    '''
    cur_dir = Path.cwd()
    source_dir = cur_dir / 'src' / 'zlib-ng'
    build_dir = cur_dir / 'zlib_build'
    build_dir.mkdir(exist_ok=True)
    os.chdir(build_dir)
    
    cmd = ['cmake', '-S', source_dir, '-B', build_dir,
        '-DZLIB_COMPAT=ON',
        '-DZLIB_ENABLE_TESTS=OFF',
        '-DBUILD_SHARED_LIBS=OFF',
        '-DCMAKE_C_FLAGS="-fPIC"',
    ]
    subprocess.run(cmd)
    subprocess.run(['cmake', '--build', build_dir, '-v', '--config', 'Release'])
    os.chdir(cur_dir)
    
    objs = [str(build_dir / 'libz.a')]
    if sys.platform == 'win32':
        objs = [str(build_dir / 'Release' / 'zlibstatic.lib'),
                ]
    
    return str(build_dir), objs

include_dir, zlib  = build_zlib()

ext = cythonize([
    Extension('pybcf.reader',
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        sources=['src/pybcf/reader.pyx',
            'src/gzstream.cpp',
            'src/bcf.cpp',
            'src/index.cpp',
            'src/header.cpp',
            'src/info.cpp',
            'src/sample_data.cpp',
            'src/variant.cpp'],
        extra_objects=zlib,
        include_dirs=['src', include_dir],
        language='c++'),
    ])

setup(name='pybcf',
    description='Package for loading data from bcf files',
    long_description=io.open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    version='1.3.6',
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