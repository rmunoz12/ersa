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
    Class StreamJSON
    ----------------
    Acts as an iterable over a function f, specifically to
    stream a list of JSON objects to json.dump() or json.dumps().

    Parameters
    ----------
    f : (**kwargs) -> generator[T]
        function that yields T

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
