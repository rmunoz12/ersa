"""Plot Number of Meioses Histogram"""
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license

import matplotlib.pyplot as plt
from sys import argv

def _read_outfile(path):
    """
    Reads an output file at path and extracts d for each related pair, and the number of unrelated pairs

    NOTE: Currently, this parses the output format of the original ERSA code, in case we decide to change our
    file format later.
    """
    d_list = []
    no_rel_count = 0
    outfile = open(path)
    lines = [[val for val in line.split()] for line in outfile]
    outfile.close()
    for line in lines:
        if line[0].startswith(("#","individual_1")):
            continue
        d = line[3]
        if d == "no_sig_rel":
            no_rel_count += 1
        else:
            d_list.append(int(d))
    return d_list, no_rel_count

def plot_histogram(d_list, no_rel_count):
    """
    Plots the optimal d for related pairs, and the number of related and unrelated pairs.
    """
    fig = plt.figure()
    dhist = fig.add_subplot(121)
    dhist.hist(d_list, bins = range(0, max(d_list) + 2), align = "left")
    dhist.set_title("Predicted Number of Meioses")
    dhist.set_xlabel("Number of Meioses")
    dhist.set_xticks(range(0, max(d_list) + 2))
    dhist.set_ylabel("Frequency")

    dfreq = fig.add_subplot(122)
    dfreq.bar([1, 2], [len(d_list), no_rel_count], width = 1.0, align = "center")
    dfreq.set_title("Significant Pairs")
    dfreq.set_xticks([1, 2])
    dfreq.set_xticklabels(("Rel", "No_Rel"))
    plt.tight_layout()
    plt.show()

def main():
    argc = len(argv)
    assert 1 < argc < 3
    path = argv[1]
    d_list, no_rel_count = _read_outfile(path)
    plot_histogram(d_list, no_rel_count)

if __name__ == '__main__':
    main()