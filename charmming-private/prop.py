#!/usr/bin/env python

import sys
from pychm.scripts.getprop import getProp

inpFileName = sys.argv[1]
outFileName = sys.argv[2]
props = sys.argv[3:] + ['avertime']

inpFile = open(inpFileName)
taco = getProp(inpFile,*props)
inpFile.close()

outLabels = []
outData = []

outLabels.append('avertime')
outData.append(map(lambda x: '%10.5f' % x,taco['avertime']))
for key in taco.keys():
    if key.endswith('time'):
        continue
    outLabels.append(key)
    outData.append(map(lambda x: '%10.5f' % x,taco[key]))

outData = map(None,*outData)

rawOutput = []
rawOutput.append('    '.join(outLabels))
for line in outData:
    rawOutput.append('    '.join(line))
rawOutput = '\n'.join(rawOutput)

outFile = open(outFileName,'w')
outFile.write(rawOutput)
outFile.close()
