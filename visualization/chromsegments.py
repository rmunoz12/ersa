import sys, os, random
import matplotlib.pyplot as plt
import numpy as np
from ersa.parser import read_matchfile, SharedSegment

def get_chrom_lengths(path):
    filename = os.path.join(os.path.dirname(__file__), path)
    chrom_dict = {}
    with open(filename) as chromfile:
        for line in [[val for val in line.split()] for line in chromfile]:
            chrom_dict[int(line[0])] = (float(line[1]), float(line[2]))
    return chrom_dict

def get_pair_dict(path):
    """
    pair_dict = {indiv1:indiv2: [(chrom, start, stop]}

    *no filtering
    *show both background and IBD matches
    """
    s_list = read_matchfile(path)
    pair_dict = {}
    for seg in s_list:
        assert isinstance(seg, SharedSegment)
        assert seg.lengthUnit == "cM"  # TODO add conversion to cM instead of assertion

        pair = seg.indivID1 + ":" + seg.indivID2
        if pair_dict.get(pair):
            pair_dict[pair].append((seg.chrom, seg.bpStart, seg.bpEnd))
        else:
            pair_dict[pair] = [(seg.chrom, seg.bpStart, seg.bpEnd)]

    return pair_dict

chrom_size = 10
chrom_sep = 5
num_autosomes = 22
def create_chrom_plot(pair_dict, chrom_dict, pair = "unspecified"):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # draw blank chromosomes
    for chrom in range(1, num_autosomes + 1):
        ax.broken_barh([(0, chrom_dict[chrom][0])], (-1 * (chrom_sep + (chrom_size + chrom_sep) * (chrom - 1)), -chrom_size), facecolors='white')

    # select pair if unspecified
    if pair == "unspecified":
        pair = random.choice(list(pair_dict.keys()))
    elif pair not in pair_dict:
        pair = pair.split(':')[1] + ':' + pair.split(':')[0]
        print(pair)
        if pair not in pair_dict:
            print("No matching entry")
            return

    # draw matching segments
    for chrom in range(1, num_autosomes + 1):
        segs = []
        for seg in [i for i in pair_dict[pair] if i[0] == chrom]:
            total_len = chrom_dict[chrom][0]
            print("chrom " + str(chrom) + " (" + str(total_len) + "):", seg[1], "-", seg[2], "|", seg[1]/total_len, "-", seg[2]/total_len)
            segs.append((seg[1], seg[2] - seg[1])) #start, width
        ax.broken_barh(segs, (-1 * (chrom_sep + (chrom_size + chrom_sep) * (chrom - 1)), -chrom_size), facecolors='blue')

    # limit chrom borders (y-axis)
    ax.set_ylim(-1 * (chrom_sep + 22 * (chrom_sep + chrom_size)), 0)
    ax.set_yticks(np.arange(-1 * (chrom_sep + chrom_size / 2), -1 * num_autosomes * (chrom_sep + chrom_size), -1 * (chrom_sep + chrom_size)))
    ax.set_yticklabels(range(1, 23))
    ax.set_ylabel('Chromosome Number')

    ax.set_xticks(np.arange(0, 3e8, 4e7))
    ax.set_xticklabels(np.arange(0, 3e2, 4e1))
    ax.set_xlabel('Chromosome Length (Mbp)')

    plt.show()

def main():
    args = sys.argv

    h = 10  #high thres (cM)
    t = 2.5
    pair_dict = get_pair_dict(args[1])
    chrom_dict = get_chrom_lengths('autosomelengths.txt')
    create_chrom_plot(pair_dict, chrom_dict)

if __name__ == '__main__':
    main()