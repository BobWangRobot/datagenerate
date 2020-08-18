from distutils.core import setup, Extension
setup(name='libsvm',
      version='3.17',
      package_dir = {'' : 'python'},
      py_modules=['svm', 'svmutil'],
      ext_modules=[Extension('libsvm', ['svm.cpp'])],
      )
