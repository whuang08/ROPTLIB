from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

print("For this compilation to work the following line has to be commented out; #define abs(x) ((x) >= 0 ? (x) : -(x))")
ext_modules = [Extension("PySimpleExample",
                     ["PySimpleExample.pyx"],
                     language='c++',
                         extra_compile_args=["-std=c++11", "-fPIC"],
                         extra_link_args=["-std=c++11"],
                    extra_objects=["libropt.so"], include_dirs=["/home/lucio/mypyprojects/ROPTLIB_2017-02-26",
                                                                     "/home/lucio/mypyprojects/ROPTLIB_2017-02-26/cwrapper/blas",
                                                                     "/home/lucio/mypyprojects/ROPTLIB_2017-02-26/cwrapper/lapack",
                                                                     "/home/lucio/mypyprojects/ROPTLIB_2017-02-26/Manifolds",
                                                                     "/home/lucio/mypyprojects/ROPTLIB_2017-02-26/Manifolds/Stiefel",
                                                                     "/home/lucio/mypyprojects/ROPTLIB_2017-02-26/Manifolds/Euclidean",
                                                                     "/home/lucio/mypyprojects/ROPTLIB_2017-02-26/Problems",
"/home/lucio/mypyprojects/ROPTLIB_2017-02-26/Solvers",
                                                                     "/home/lucio/mypyprojects/ROPTLIB_2017-02-26/Problems/StieBrockett",
                                                                     "/home/lucio/mypyprojects/ROPTLIB_2017-02-26/Others"]
                     )]

#setup(ext_modules=cythonize("PyManifoldPolLDR.pyx"))
setup(
  name = 'PySimpleExample',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)