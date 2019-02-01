from distutils.core import setup, Extension
import os



setup(name='websky_model',
      version='0.1',
      description='Data Model for the WebSky sims',
      url='https://github.com/simonsobs/websky_model/',
      license='BSD-2-Clause',
      packages=['websky_model'],
      package_dir={'websky_model':'websky_model'},
      zip_safe=False)
