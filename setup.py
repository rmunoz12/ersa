from setuptools import setup, find_packages

DESCRIPTION = "Implementation of Huff et. al. (2011) Estimation of Recent Shared Ancestry"
LONG_DESCRIPTION = DESCRIPTION
NAME = "ersa"
AUTHOR = "Richard Munoz, Jie Yuan, Yaniv Erlich"
AUTHOR_EMAIL = "rmunoz@columbia.edu, jyuan@columbia.edu, yaniv@cs.columbia.edu"
MAINTAINER = "Richard Munoz"
MAINTAINER_EMAIL = "rmunoz@columbia.edu"
DOWNLOAD_URL = 'http://github.com/rmunoz12/ersa'
LICENSE = 'GNU GPL v3'

VERSION = '1.0.0'

setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      url=DOWNLOAD_URL,
      download_url=DOWNLOAD_URL,
      license=LICENSE,
      packages=find_packages(),
      test_suite='tests',
      entry_points = {
          "console_scripts": ['ersa = ersa.ersa:main']
      },
      install_requires=['sqlalchemy', 'inflect', 'pytest', 'scipy'],
      classifiers=['Development Status :: 4 - Beta',
                   'Programming Language :: Python :: 3.4',
                   'License :: OSI Approved :: GNU GPL v3',
                   'Operating System :: OS Independent',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: Bio-Informatics']
      )
