#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt
import re
import sys
import itertools
import pickle

ap = argparse.ArgumentParser(description="Plot throughput and cumulative sequence by time.")

ap.add_argument("sam", help="Input sam file. This must be generated with blasr supplied by https://github.com/mchaisso/blasr, and ran with the options '-samqv IPD -sam -clipping soft -bestn 1 ' in addition to the bas.h5 file and reference.")

ap.add_argument("output", help="Base name for output file.")
ap.add_argument("--frameRate", help="Frame rate (75).", default=75, type=int)
ap.add_argument("--timeBin", help="Bin into this number of seconds.", default=60, type=int)
ap.add_argument("--csv", help="Write values to a CSV", default=None)
ap.add_argument("--pickle", help="Write files to a pickle.", default=None)
args = ap.parse_args()
samFile = open(args.sam)
alignTimes = []
index = 0
alnTotal = 0
endTimes = []
seqLengths = []
for line in samFile:
    if (line[0] == '@'):
        continue
    else:
        vals = line.split()
        foundIPD = False
        cigar = vals[5]
        ops = re.split("\d+", cigar)[1:]
        lengths = [int(i) for i in re.split("[MIDNSHP=X]", cigar)[:-1]]
        softClipPrefix = 0
        seqLen = len(vals[9])
        softClipSuffix = seqLen

        firstM = 0
        while (firstM < len(lengths) and ops[firstM] != 'M'):
            firstM+=1
        lastM = len(lengths) -1

        while (lastM > firstM and ops[lastM] != 'M'):
            lastM -= 1
            
        if (ops[firstM-1] == 'S'):
            softClipPrefix = lengths[firstM-1]
        if (lastM < len(lengths)-1 and ops[lastM+1] == 'S'):
            softClipSuffix = seqLen - lengths[lastM+1]

        alnTotal += softClipSuffix - softClipPrefix
        for i in range(9,len(vals)):
            if ("ip:Z:S" in vals[i]):
                foundIPD = True
                break
        if (foundIPD == True):
            ipdStr=vals[i][6:]
            ipd = [int(v) for v in ipdStr.split(",")]
            startTime = np.cumsum(np.asarray(ipd))
            startTime = startTime / args.frameRate
            alignTime = startTime[softClipPrefix:softClipSuffix] + startTime[softClipPrefix]
            alignTimes.append(alignTime)
            endTimes.append(alignTime[-1])
            seqLengths.append(softClipSuffix - softClipPrefix)
        else:
            print "no ipd for " + str(softClipSuffix - softClipPrefix) + " " + vals[0]
        index += 1
        if (index % 1000 == 0):
            sys.stderr.write("processed {}\n".format(index))

        
maxTime = 0
for i in range(0,len(alignTimes)):
    if (len(alignTimes[i]) > 0):
        maxTime = max(alignTimes[i].max(), maxTime)

print "max time is: " + str(maxTime)
print "nbins: " + str(maxTime / args.timeBin)
print "align total: " + str(alnTotal)
ninc = np.zeros(maxTime/ args.timeBin + 1)
total = 0
for i in range(0,len(alignTimes)):
#    print "len align times /attotal"
#    print len(alignTimes[i])
    atTotal = 0
    for time,count in itertools.groupby(alignTimes[i]/args.timeBin):
        l = list(count)
#        print "{} : {} ({})".format(time, l[0], len(l))
#        print str(l)
        n = len(l)
        total += n
        atTotal += n
        ninc[time] += n
#    print atTotal


print "total " + str(total)
print len(alignTimes)
inccs = np.cumsum(ninc)
print ninc
print inccs
if (args.pickle is not None):
    pf = open(args.pickle, 'w')
    pickle.dump(ninc, pf)
    pickle.dump(inccs, pf)
    pickle.dump(endTimes, pf)
    pickle.dump(seqLengths, pf)
    pf.close()


plt.plot(ninc)
plt.xlabel("Time (per {} sec)".format(args.timeBin))
plt.ylabel("Bin # bases")
plt.title("Incorporation rate.")
plt.savefig(args.output + ".rate.png")
plt.clf()
plt.plot(inccs)
plt.xlabel("Time (per {} sec)".format(args.timeBin))
plt.ylabel("Total # bases")
plt.title("Incorporation rate.")
plt.savefig(args.output + ".total.png")
plt.clf()
plt.plot(np.asarray(endTimes)/args.timeBin, np.asarray(seqLengths), 'bo')
plt.xlabel("Time (per {} sec)".format(args.timeBin))
plt.ylabel("Seq length (bases)")
plt.title("End times for sequence length")
plt.savefig(args.output + ".lengths.png")

if (args.csv is not None):
    totals = np.ndarray([2,len(ninc)])
    totals[0] = ninc
    totals[1] = inccs
    np.savetxt(args.csv, totals, fmt="%d", delimiter=",")
    
