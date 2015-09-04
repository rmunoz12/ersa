"""Reading and parsing *.match file from Germline"""
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license

class SharedSegment:
    """
    Class SharedSegment
    -------------------

    Stores named values for a matchfile line, see "Output" from
    http://www1.cs.columbia.edu/~gusev/germline/

    Assumes that each value in the inputed param_list is a string.
    """
    def __init__(self, param_list):
        self.familyID1 = int(param_list[0])
        self.indivID1 = param_list[1]
        self.familyID2 = int(param_list[2])
        self.indivID2 = param_list[3]
        self.chrom = int(param_list[4])
        self.bpStart = int(param_list[5])
        self.bpEnd = int(param_list[6])
        self.snpStart = param_list[7]
        self.snpEnd = param_list[8]
        self.totalSNP = int(param_list[9])
        self.length = float(param_list[10])
        self.lengthUnit = param_list[11]
        self.mismatchSNP = int(param_list[12])
        self.ind1homozygous = int(param_list[13])
        self.ind2homozygous = int(param_list[14])


def read_matchfile(path):
    """
    Reads a match file at path and returns a list of SharedSegments.
    """
    s_list = []
    matchfile = open(path)
    lines = [[val for val in line.split()] for line in matchfile]
    matchfile.close()
    for line in lines:
        segment = SharedSegment(line)
        s_list.append(segment)
    return s_list


def get_pair_dict(path):
    """
    Reads and collapses the input data into a dictionary with entries:

    {pair_id: (n, s)}

    pair_id = unique id for each pair of individuals compared
    n = number of shared segments
    s = list of shared segment lengths (in cM)
    """
    s_list = read_matchfile(path)
    pair_dict = {}
    for seg in s_list:
        assert isinstance(seg, SharedSegment)
        assert seg.lengthUnit == "cM"  # TODO add conversion to cM instead of assertion
        pair_id = seg.indivID1
        if seg.indivID1 < seg.indivID2:
            pair_id += ":" + seg.indivID2
        else:
            pair_id = seg.indivID2 + ":" + pair_id
        if pair_dict.get(pair_id):
            n = pair_dict[pair_id][0] + 1
            shared_segs = pair_dict[pair_id][1]
            shared_segs.append(seg.length)
            pair_dict[pair_id] = (n, shared_segs)
        else:
            pair_dict[pair_id] = (1, [seg.length])
    return pair_dict


def main():
    paired_data = get_pair_dict("../test_data/test_parser.match")
    for k, v in paired_data.items():
        print("pair_id: {}\t(n, s): {}".format(k, v))


if __name__ == '__main__':
    main()
