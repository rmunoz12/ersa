"""
Genomic Region Masking
"""
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license


"""
b : int
    regions must extend this far beyond each end of a masked region
"""
b = 1 * 10 ** 6


"""
masked : list[list[Tuple[int, int]]
    genomic regions that should be masked, see Huff et al. 2014 (Table 3);
    in total there are 14 regions and 119.92 cM;
    all masked regions are of at least 5 cM
"""
masked = [[] for x in range(23)]
masked[9].append((38293483, 72605261, 8.15))
masked[8].append((10428647, 13469693, 7.96))
masked[21].append((16344186, 19375168, 6.91))
masked[10].append((44555093, 53240188, 7.58))
masked[22].append((16051881, 25095451, 20.82))
masked[2].append((85304243, 99558013, 6.53))
masked[1].append((118434520, 153401108, 9.95))
masked[15].append((20060673, 25145260, 10.46))
masked[17].append((77186666, 78417478, 5.66))
masked[15].append((27115823, 30295750, 9.29))
masked[17].append((59518083, 64970531, 6.23))
masked[2].append((132695025, 141442636, 9.16))
masked[16].append((19393068, 24031556, 6.18))
masked[2].append((192352906, 198110229, 5.04))


def mask_input_segs(segs, t):
    """
    Modifies input segment lengths that fall into
    masked regions.

    Parameters
    ----------
    segs : list[SharedSegment]
        original input SharedSegments

    t : float
        Filter out results less than t (in cM)

    Returns
    -------
    new_segs : list[SharedSegment]
        modified segs list, with modified lengths and
        fully masked segments removed
    """
    new_segs = []
    for s in segs:
        regions = masked[s.chrom]
        for r in regions:
            l = r[0]  # low
            h = r[1]  # high
            mask_cm = r[2]
            # todo check if l - b or h + b fall off a chromosome?
            if s.bpStart > l - b and s.bpEnd < h + b:
                s.length = 0
                break
            elif s.bpStart <= l - b and s.bpEnd >= h + b:
                # subtract len but keep start & end
                # for visual display
                s.length -= mask_cm
                break
            elif l <= s.bpStart <= h < s.bpEnd:
                # truncate, changing start
                # for the visual display
                ratio = (s.bpEnd - h) / (s.bpEnd - s.bpStart)
                s.length *= ratio
                s.length = round(s.length, 2)
                s.bpStart = h
                break
            elif h >= s.bpEnd >= l > s.bpStart:
                # truncate, changing end
                # for the visual display
                ratio = (l - s.bpStart) / (s.bpEnd - s.bpStart)
                s.length *= ratio
                s.length = round(s.length, 2)
                s.bpEnd = l
                break
        if s.length >= t:
            new_segs.append(s)

    return new_segs


def total_masked():
    """
    Returns the sum of the lengths of masked regions in cM

    Returns
    -------
    m : float
    """
    m = 0
    for c in masked:
        for r in c:
            m += r[2]
    return m
