from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension("GalaxyvsContaminant", ["GalaxyvsContaminant.pyx"],
	libraries=["m"],
	extra_compile_args = ["-O3", "-ffast-math", "-march=native", "-fopenmp" ,"-Wno-cpp"],
	include_dirs=[numpy.get_include()],
	extra_link_args=['-fopenmp']),
	Extension("template_fit", ["template_fit.pyx"],
	libraries=["m"],
	extra_compile_args = ["-ffast-math","-Wno-cpp"],
	include_dirs=[numpy.get_include()]),
	Extension("rho_r", ["rho_r.pyx"],
	libraries=["m"],
	extra_compile_args = ["-ffast-math","-Wno-cpp"],
	include_dirs=[numpy.get_include()])

]
setup(
	ext_modules=cythonize(extensions),
)
