import os
import sys

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

#filename="/Users/jyuan/Desktop/generated.match"
segmentList = []
filename = str(sys.argv[1])
#print "path to file: ", filename
matchfile = open(filename, "r")

for line in matchfile.readlines():
    parameterList = line.strip('\n').replace(' ','\t').split('\t')
    segment = SharedSegment(parameterList)
    segmentList.append(segment)
matchfile.close()

for segment in segmentList:
    print segment