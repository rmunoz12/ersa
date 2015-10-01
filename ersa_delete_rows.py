#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license

from ersa.dbmanager import DbManager
from argparse import ArgumentParser

def get_args():
    p = ArgumentParser(description="Hard delete rows marked as \"deleted\" in a database created with ersa")
    p.add_argument("db", help="database to delete rows from")
    args = p.parse_args()
    return args

def hard_delete_rows():
    args = get_args()
    with DbManager(args.db) as db:
        db.delete()



if __name__ == '__main__':
    hard_delete_rows()
