import sys

"""
Class SharedSegment
-------------------

Stores named values for a matchfile line, see "Output" from
http://www1.cs.columbia.edu/~gusev/germline/

Assumes that each value in the parameter list
is a string.
"""
class SharedSegment:
    def __init__(self, parameterList):
        self.familyID1 = parameterList[0]
        self.indivID1 = parameterList[1]
        self.familyID2 = parameterList[2]
        self.indivID2 = parameterList[3]
        self.chrom = parameterList[4]
        self.bpStart = parameterList[5]
        self.bpEnd = parameterList[6]
        self.snpStart = parameterList[7]
        self.snpEnd = parameterList[8]
        self.totalSNP = parameterList[9]
        self.length = parameterList[10]
        self.lengthUnit = parameterList[11]
        self.mismatchSNP = parameterList[12]
        self.ind1homozygous = parameterList[13]
        self.ind2homozygous = parameterList[14]


def read_matchfile(path):
    sList = []
    matchfile = open(path)
    lines = [[val for val in line.split()] for line in matchfile]
    matchfile.close()
    for line in lines:
        segment = SharedSegment(line)
        sList.append(segment)
    return sList


def get_pair_dict(path):
    """
    Reads and collapses the input data into a dictionary with entries:

    {pair_id: (n, s)}

    pair_id = unique id for each pair of individuals compared
    n = number of shared segments
    s = list of shared segment lengths (in cM)
    """
    sList = read_matchfile(path)
    pair_dict = {}
    for seg in sList:
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
    paired_data = get_pair_dict("../test_data/generated.match")
    for k, v in paired_data.items():
        print("paird_id: {}\t(n, s): {}".format(k, v))


if __name__ == '__main__':
    main()
