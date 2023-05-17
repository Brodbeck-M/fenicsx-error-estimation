# Copyright 2019-2020, Jack S. Hale, Raphaël Bulle
# SPDX-License-Identifier: LGPL-3.0-or-later
import os
import sys
import platform
import subprocess

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

if sys.version_info < (3, 7):
    print("Python 3.7 or higher required, please upgrade.")
    sys.exit(1)

VERSION = "0.5.0"
URL = ""

REQUIREMENTS = [
    "numpy",
    "scipy",
    "fenics-dolfinx>=0.5.0",
]

AUTHORS = """\
Raphael Bulle
Jack S. Hale
"""

CLASSIFIERS = """\
Development Status :: 5 - Production/Stable
Intended Audience :: Developers
Intended Audience :: Science/Research
License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)
Operating System :: POSIX
Operating System :: POSIX :: Linux
Operating System :: MacOS :: MacOS X
Programming Language :: Python
Programming Language :: Python :: 3
Programming Language :: Python :: 3.7
Programming Language :: Python :: 3.8
Topic :: Scientific/Engineering :: Mathematics
"""


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the"
                               + "following extensions:"
                               + ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            raise RuntimeError("Windows not supported")
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(
            self.get_ext_fullpath(ext.name)))
        cmake_args = [f'-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}',
                      f'-DPYTHON_EXECUTABLE={sys.executable}']

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += [
                f'-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}']
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j3']

        env = os.environ.copy()
        import pybind11
        env['pybind11_DIR'] = pybind11.get_cmake_dir()
        env['CXXFLAGS'] = f'{env.get("CXXFLAGS", "")} -DVERSION_INFO=\\"{self.distribution.get_version()}\\"'

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args,
                              cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args,
                              cwd=self.build_temp, env=env)

setup(
    name="fenicsx_error_estimation",
    description="Implicit and a posteriori error estimates in FEniCSx",
    version=VERSION,
    author=AUTHORS,
    classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
    license="LGPL version 3 or later",
    author_email="mail@jackhale.co.uk",
    maintainer_email="mail@jackhale.co.uk",
    url=URL,
    packages=["fenicsx_error_estimation"],
    ext_modules=[CMakeExtension('fenicsx_error_estimation.cpp')],
    cmdclass=dict(build_ext=CMakeBuild),
    platforms=["Linux", "Mac OS-X", "Unix"],
    install_requires=REQUIREMENTS,
    zip_safe=False
)