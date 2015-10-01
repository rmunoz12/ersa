#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license

from ersa.dbmanager import DbManager
from argparse import ArgumentParser
import sys


def get_args():
    p = ArgumentParser(description="Hard delete rows marked as \"deleted\" in a database created with ersa")
    p.add_argument("db", help="database to delete rows from")
    args = p.parse_args()
    return args


def query_yes_no(question, default="yes"):
    """
    Ask a yes/no question via raw_input() and return their answer.

    question : str
        prompt presented to the user

    default : str | None
        The presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".

    References
    ----------
    http://stackoverflow.com/questions/3041986/python-command-line-yes-no-input
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


def hard_delete_rows():
    args = get_args()

    question = "Preparing to delete rows marked as deleted from \"" + \
               args.db + "\". Proceed?"

    if query_yes_no(question, default="no"):
        with DbManager(args.db) as db:
            db.delete()
    else:
        print("No rows deleted. Exiting..")


if __name__ == '__main__':
    hard_delete_rows()
