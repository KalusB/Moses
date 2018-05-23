from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

GalaxyvsContaminant = [
    Extension("GalaxyvsContaminant", ["GalaxyvsContaminant.pyx"],
	libraries=["m"],
	extra_compile_args = ["-ffast-math","-Wno-cpp"],
    include_dirs=[numpy.get_include()])
]
setup(
	ext_modules=cythonize(GalaxyvsContaminant),
)
