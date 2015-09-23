# ersa
`ersa` estimates the combined number of generations between pairs of individuals.  Takes a [Germline](http://www1.cs.columbia.edu/~gusev/germline/) matchfile as input.

`ersa` is an implementation of [Huff et. al. (2011) Maximum-Likelihood estimation of recent shared ancenstry (ERSA)](http://genome.cshlp.org/content/21/5/768.full).

## Example
Requires a germline matchfile as input.  By default, summary results are sent to stdout, but can also be directed to a file.  Alternatively, full results can be directed to a database, as in the example below.

    $ python3 ersa.py -D "sqlite:///ersa_results.db" input.match

This creates `ersa_results.db` in current directory.

For additional options, use

    $ python3 ersa.py -h
