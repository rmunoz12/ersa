from setuptools import setup, find_packages

DESCRIPTION = "Implementation of Huff et. al. (2011) Estimation of Recent Shared Ancestry "
LONG_DESCRIPTION = "`ersa` estimates the combined number of generations between pairs of " \
                   "individuals using a " \
                   "`Germline <http://www1.cs.columbia.edu/~gusev/germline/>`_ " \
                   "matchfile as input.  It is an implementation of " \
                   "`Huff et. al. (2011) Maximum-Likelihood estimation of recent shared ancenstry (ERSA) <http://genome.cshlp.org/content/21/5/768.full>`_ " \
                   "and `Li et. al. (2014) Relationship Estimation from Whole-Genome Sequence Data <http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004144>`_ ."
NAME = "ersa"
AUTHOR = "Richard Munoz, Jie Yuan, Yaniv Erlich"
AUTHOR_EMAIL = "rmunoz@columbia.edu, jyuan@columbia.edu, yaniv@cs.columbia.edu"
MAINTAINER = "Richard Munoz"
MAINTAINER_EMAIL = "rmunoz@columbia.edu"
DOWNLOAD_URL = 'http://github.com/rmunoz12/ersa'
LICENSE = 'GNU General Public License v3 (GPLv3)'

VERSION = '1.1.1'

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
                   'Intended Audience :: Science/Research',
                   'Programming Language :: Python :: 3.4',
                   'Programming Language :: Python :: 3 :: Only',
                   'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                   'Operating System :: OS Independent',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Topic :: Scientific/Engineering :: Information Analysis',
                   'Topic :: Sociology :: Genealogy']
      )
