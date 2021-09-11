from distutils.core import setup
from Cython.Build import cythonize

ext_options = {"compiler_directives": {"profile": True}, "annotate": True}

setup(ext_modules = cythonize("coherence_gain_t_finite_detailed.pyx", **ext_options)
)
