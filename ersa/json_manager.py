"""
JSON Output Handling
"""
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license


class StreamJSON(list):
    """
    Parameters
    ----------
    f
        function that yields a generator and takes kwwargs

    Notes
    -----
    Converts f into an iterable higher order function

    References
    ----------
    http://stackoverflow.com/questions/21663800/python-make-a-list-generator-json-serializable
    """
    def __init__(self, f, **kwargs):
        super().__init__()
        self.gen = self._increase_order(f, **kwargs)

    def __iter__(self):
        return self.gen()

    def __len__(self):
        return 1

    @staticmethod
    def _increase_order(f, **kwargs):
        def helper():
            return f(**kwargs)
        return helper
