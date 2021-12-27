from sys import argv
import uc
from os import exit
from fasta import GetSizeFromLabel

FileName = argv[1]


def GetSampleId(Label):
    Fields = Label.split(";")
    for Field in Fields:
        if Field.startswith("barcodelabel="):
            return Field[13:]
    print("barcodelabel= not found in read label " + str(Label))
    exit(1)


def OnRec():
    global OTUs, Samples, OTUTable
    if uc.Type != 'H':
        return

    OTUId = uc.TargetLabel
    if OTUId not in OTUIds:
        OTUIds.append(OTUId)
        OTUTable[OTUId] = {}

    SampleId = GetSampleId(uc.QueryLabel)
    if SampleId not in SampleIds:
        SampleIds.append(SampleId)

    N = GetSizeFromLabel(uc.QueryLabel, 1)
    try:
        OTUTable[OTUId][SampleId] += N
    except BaseException:
        OTUTable[OTUId][SampleId] = N


OTUIds = []
SampleIds = []
OTUTable = {}

uc.ReadRecs(FileName, OnRec)

s = "OTUId"
for SampleId in SampleIds:
    s += "\t" + SampleId
print(s)

for OTUId in OTUIds:
    s = OTUId
    for SampleId in SampleIds:
        try:
            n = OTUTable[OTUId][SampleId]
        except BaseException:
            n = 0
        s += "\t" + str(n)
    print(s)
