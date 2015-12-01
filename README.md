# ersa
`ersa` estimates the combined number of generations between pairs of individuals using a [Germline](http://www1.cs.columbia.edu/~gusev/germline/) matchfile as input.  It is an implementation of [Huff et. al. (2011) Maximum-Likelihood estimation of recent shared ancenstry (ERSA)](http://genome.cshlp.org/content/21/5/768.full) and [Li et. al. (2014) Relationship Estimation from Whole-Genome Sequence Data](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004144).

Copyright (c) 2015 by
- Richard Muñoz (rmunoz@nygenome.org)
- Jie Yuan (jyuan@nygenome.org)
- Yaniv Erlich (yaniv@nygenome.org)

License: GNU GPL v3 (see LICENSE.txt)

## Install
First, install Python 3.4 or greater (and setuptools, if necessary).  Then clone from github:

    git clone https://github.com/rmunoz12/ersa.git
    cd ersa
    sudo python3 setup.py install

This will add `ersa` to your path and download additional required python packages.

## Example
Requires a Germline matchfile as input.  By default, results are sent to standard out in JSON format, but can also be directed to a file.  Alternatively, results can be directed to a database, as in the example below.

    $ ersa -D "sqlite:///ersa_results.db" input.match

This creates `ersa_results.db` in current directory.

For additional options, use

    $ ersa --help

## Notes
On inserting results into a database, if a comparison between a pair of individuals exists, `ersa` will mark the old result as deleted (i.e., soft delete the result).  To physically delete these old results from the database, a utlity `ersa_delete_rows` is also provided:

    $ ersa_delete_rows "sqlite:///ersa_results.db"

