# ersa
`ersa` estimates the combined number of generations between pairs of individuals using a [Germline](http://www1.cs.columbia.edu/~gusev/germline/) matchfile as input.  It is an implementation of [Huff et. al. (2011) Maximum-Likelihood estimation of recent shared ancenstry (ERSA)](http://genome.cshlp.org/content/21/5/768.full).

Copyright (c) 2015 by
- Richard Muñoz (rmunoz@nygenome.org)
- Jie Yuan (jyuan@nygenome.org)
- Yaniv Erlich (yaniv@nygenome.org)

License: GNU GPL v3 (see LICENSE.txt)

## Requirements
`ersa` requires Python 3.4 or greater.  In addition, the following packages for python3 must be installed prior to the installation process:

- `numpy`
- `setuptools`

## Install
First, install python3, setuptools, and numpy.  Then clone from github:

    git clone https://github.com/rmunoz12/ersa.git
    cd ersa
    sudo python3 setup.py install

This will add `ersa` to your path and download additional required python packages.

## Example
Requires a Germline matchfile as input.  By default, summary results are sent to stdout, but can also be directed to a file.  Alternatively, full results can be directed to a database, as in the example below.

    $ ersa -D "sqlite:///ersa_results.db" input.match

This creates `ersa_results.db` in current directory.

For additional options, use

    $ ersa -h
