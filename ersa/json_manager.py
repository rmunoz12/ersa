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
    def __init__(self, gen, **kwargs):
        super().__init__()
        self.gen = self._increase_order(gen, **kwargs)

    def __iter__(self):
        return self.gen()

    def __len__(self):
        return 1

    @staticmethod
    def _increase_order(f, **kwargs):
        def helper():
            return f(**kwargs)
        return helper
