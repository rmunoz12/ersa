""" Reading and parsing *.match file from Germline """
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license

from ersa.mask import mask_input_segs


class SharedSegment:
    """
    Structure that stores named values for a matchfile line, see
    "Output" from http://www1.cs.columbia.edu/~gusev/germline/

    Parameters
    ----------
    param_list : list[str]
        list of strings in the order of a germline matchfile.

    Notes
    -----
    SharedSegments are ordered by the length parameter.
    """
    def __init__(self, param_list):
        assert type(param_list) == list
        assert len(param_list) == 15
        # self.familyID1 = int(param_list[0])
        self.indivID1 = param_list[1]
        # self.familyID2 = int(param_list[2])
        self.indivID2 = param_list[3]
        self.chrom = int(param_list[4])
        self.bpStart = int(param_list[5])
        self.bpEnd = int(param_list[6])
        # self.snpStart = param_list[7]
        # self.snpEnd = param_list[8]
        # self.totalSNP = int(param_list[9])
        self.length = float(param_list[10])
        self.lengthUnit = param_list[11]
        # self.mismatchSNP = int(param_list[12])
        # self.ind1homozygous = int(param_list[13])
        # self.ind2homozygous = int(param_list[14])

    def __lt__(self, other):
        return self.length < other.length


def read_matchfile(path, haploscores=False):
    """
    Reads a matchfile at path and yields SharedSegments.

    Parameters
    ----------
    path : str

    haploscores : bool
        True if the input matchfile contains haploscores in an
        extra column at the end of each line. These scores
        are discarded.

    Returns
    -------
    segment : generator[SharedSegment]
    """
    with open(path) as matchfile:
        for line in matchfile:
            split_line = [val for val in line.split()]
            if haploscores:
                del split_line[-1:]
            segment = SharedSegment(split_line)
            yield segment


def get_pair_dict(path, t, user=None, haploscores=False, nomask=False):
    """
    Reads from path and collapses the input data into a dictionary
    mapping pairs to SharedSegments.

    Parameters
    ----------
    path : str

    t : float
        Filter out results less than t (in cM)

    user : str | None
        filter input by user identification

    haploscores : bool
        True if the input matchfile contains haploscores in an
        extra column at the end of each line.

    Returns
    -------
    pair_dict: dict[str: list[SharedSegments]]
        Each list of SharedSegments is sorted for processing by
        ersa_LL.estimate_relation()
    """
    s_list = read_matchfile(path, haploscores)
    pair_dict = {}
    for seg in s_list:
        assert isinstance(seg, SharedSegment)
        assert seg.lengthUnit == "cM"  # TODO add conversion to cM instead of assertion
        if seg.length < t:  # Note: seg.length > h filtered only for background parameters
            continue
        if user and seg.indivID1 != user and seg.indivID2 != user:
            continue

        pair_id = seg.indivID1
        if seg.indivID1 < seg.indivID2:
            pair_id += ":" + seg.indivID2
        else:
            pair_id = seg.indivID2 + ":" + pair_id
        if pair_dict.get(pair_id):
            pair_dict[pair_id].append(seg)
        else:
            pair_dict[pair_id] = [seg]

    for pair, segs in pair_dict.items():
        if not nomask:
            segs = mask_input_segs(segs, t)
            pair_dict[pair] = segs
        segs.sort()

    return pair_dict
