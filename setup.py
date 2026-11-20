
from pathlib import Path
import platform
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import subprocess
import sys

from Cython.Build import cythonize

ROOT = Path(__file__).resolve().parent

EXTRA_COMPILE_ARGS = []
if sys.platform == 'linux':
    EXTRA_COMPILE_ARGS += ['-std=c++11', '-O2']
elif sys.platform == "darwin":
    EXTRA_COMPILE_ARGS += ["-stdlib=libc++", "-std=c++11", "-O2"]
elif sys.platform == "win32":
    EXTRA_COMPILE_ARGS += ['/std:c++14', '/O2']

if platform.machine() == 'x86_64' and sys.platform == "linux":
    EXTRA_COMPILE_ARGS += ['-mavx', '-mavx2']

if platform.machine() == 'aarch64' and sys.platform == "linux":
    EXTRA_COMPILE_ARGS += ['-flax-vector-conversions']

def build_zlib(build_dir):
    ''' compile zlib-ng into a static library to link with
    
    Args:
        build_dir: folder to build in, which also gets the generated zlib.h
    
    Returns:
        list holding the path to the static library
    '''
    source_dir = ROOT / 'src' / 'zlib-ng'
    cmd = ['cmake', '-S', source_dir, '-B', build_dir,
        '-DZLIB_COMPAT=ON',
        '-DZLIB_ENABLE_TESTS=OFF',
        '-DBUILD_SHARED_LIBS=OFF',
        '-DCMAKE_POSITION_INDEPENDENT_CODE=ON',
    ]
    subprocess.run(cmd, check=True)
    subprocess.run(['cmake', '--build', build_dir, '-v', '--config', 'Release'],
                   check=True)
    
    if sys.platform == 'win32':
        return [str(build_dir / 'Release' / 'zlibstatic.lib')]
    return [str(build_dir / 'libz.a')]

class build_ext_subclass(build_ext):
    ''' link the vendored zlib-ng into the extension module
    
    zlib-ng is only needed to link the extension module, so compiling it here
    (rather than at import) keeps metadata-only commands such as egg_info and
    sdist from invoking cmake.
    '''
    def build_extensions(self):
        # built alongside the other temporary files, rather than in the source tree
        zlib_dir = Path(self.build_temp) / 'zlib-ng'
        zlib = build_zlib(zlib_dir)
        for ext in self.extensions:
            ext.include_dirs.append(str(zlib_dir))
            ext.extra_objects = zlib
        super().build_extensions()

extensions = [
    Extension('pybcf.reader',
        extra_compile_args=EXTRA_COMPILE_ARGS,
        sources=['src/pybcf/reader.pyx',
            'src/gzreader.cpp',
            'src/bcf.cpp',
            'src/index.cpp',
            'src/header.cpp',
            'src/info.cpp',
            'src/sample_data.cpp',
            'src/variant.cpp'],
        include_dirs=['src'],
        language='c++'),
    ]

setup(
    package_dir={'': 'src'},
    ext_modules=cythonize(extensions),
    cmdclass={'build_ext': build_ext_subclass},
    )
