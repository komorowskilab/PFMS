#
#  commoncode.py
#  ENRAGE
#

import tempfile
import shutil
import os
from os import environ
import string
import sqlite3 as sqlite
from time import strftime
from array import array
from collections import defaultdict

commoncodeVersion = 5.5
currentRDSversion = 1.1

if environ.get("CISTEMATIC_TEMP"):
    cisTemp = environ.get("CISTEMATIC_TEMP")
else:
    cisTemp = "/tmp"

tempfile.tempdir = cisTemp


def getReverseComplement(base):
    revComp = {"A": "T",
               "T": "A",
               "G": "C",
               "C": "G",
               "N": "N"
        }

    return revComp(base)


def countDuplicatesInList(listToCheck):
    tally = defaultdict(int)
    for item in listToCheck:
        tally[item] += 1

    return tally.items()


def writeLog(logFile, messenger, message):
    """ create a log file to write a message from a messenger or append to an existing file.
    """
    try:
        logfile = open(logFile)
    except IOError:
        logfile = open(logFile, "w")
    else:
        logfile = open(logFile, "a")

    logfile.writelines("%s: [%s] %s\n" % (strftime("%Y-%m-%d %H:%M:%S"), messenger, message))
    logfile.close()


def getMergedRegions(regionfilename, maxDist=1000, minHits=0, verbose=False, keepLabel=False,
                     fullChrom = False, chromField=1, scoreField=4, pad=0, compact=False,
                     doMerge=True, keepPeak=False, returnTop=0):
    """ returns a list of merged overlapping regions; 
    can optionally filter regions that have a scoreField fewer than minHits.
    Can also optionally return the label of each region, as well as the
    peak, if supplied (peakPos and peakHeight should be the last 2 fields).
    Can return the top regions based on score if higher than minHits.
    """
    regions = {}
    infile = open(regionfilename)
    lines = infile.readlines()
    hasPvalue = 0
    hasShift = 0
    if 0 < returnTop < len(lines):
        scores = []
        for line in lines:
            if line[0] == "#":
                if "pvalue" in line:
                    hasPvalue = 1

                if "readShift" in line:
                    hasShift = 1

                continue

            fields = line.strip().split("\t")
            hits = float(fields[scoreField].strip())
            scores.append(hits)

        scores.sort()
        returnTop = -1 * returnTop 
        minScore = scores[returnTop]
        if minScore > minHits:
            minHits = minScore

    mergeCount = 0
    chromField = int(chromField)
    count = 0
    #TODO: Current algorithm processes input file line by line and compares with prior lines.  Problem is it
    #      exits at the first merge.  This is not a problem when the input is sorted by start position, but in
    #      the case of 3 regions ABC that are in the input file as ACB as it goes now when processing C there
    #      will be no merge with A as B is needed to bridge the two.  When it comes time to process B it will
    #      be merged with A but that will exit the loop and the merge with C will be missed.
    for line in lines:
        if line[0] == "#":
            if "pvalue" in line:
                hasPvalue = 1

            if "readShift" in line:
                hasShift = 1

            continue

        fields = line.strip().split("\t")
        if minHits >= 0:
            try:
                hits = float(fields[scoreField].strip())
            except:
                continue

            if hits < minHits:
                continue

        if compact:
            (chrom, pos) = fields[chromField].split(":")
            if not fullChrom:
                chrom = chrom[3:]

            (front, back) = pos.split("-")
            start = int(front)
            stop = int(back)
        elif chromField > 1:
            label = string.join(fields[:chromField],"\t")
            if fullChrom:
                chrom = fields[chromField]
            else:
                chrom = fields[chromField][3:]

            start = int(fields[chromField + 1]) - pad
            stop = int(fields[chromField + 2]) + pad
        else:
            label = fields[0]
            if fullChrom:
                chrom = fields[1]
            else:
                chrom = fields[1][3:]

            start = int(fields[2]) - pad
            stop = int(fields[3]) + pad

        length = abs(stop - start)
        if keepPeak:
            peakPos = int(fields[-2 - hasPvalue - hasShift])
            peakHeight = float(fields[-1 - hasPvalue - hasShift])

        if chrom not in regions:
            regions[chrom] = []

        merged = False

        if doMerge and len(regions[chrom]) > 0:
            for index in range(len(regions[chrom])):
                if keepLabel and keepPeak:
                    (rlabel, rstart, rstop, rlen, rpeakPos, rpeakHeight) = regions[chrom][index]
                elif keepLabel:
                    (rlabel, rstart, rstop, rlen) = regions[chrom][index]
                elif keepPeak:
                    (rstart, rstop, rlen, rpeakPos, rpeakHeight) = regions[chrom][index]
                else:
                    (rstart, rstop, rlen) = regions[chrom][index]

                if regionsOverlap(start, stop, rstart, rstop) or regionsAreWithinDistance(start, stop, rstart, rstop, maxDist):
                    if start < rstart:
                        rstart = start

                    if rstop < stop:
                        rstop = stop

                    rlen = abs(rstop - rstart)
                    if keepPeak:
                        if peakHeight > rpeakHeight:
                            rpeakHeight = peakHeight
                            rpeakPos = peakPos

                    if keepLabel and keepPeak:
                        regions[chrom][index] = (label, rstart, rstop, rlen, rpeakPos, rpeakHeight)
                    elif keepLabel:
                        regions[chrom][index] = (label, rstart, rstop, rlen)
                    elif keepPeak:
                        regions[chrom][index] = (rstart, rstop, rlen, rpeakPos, rpeakHeight)
                    else:
                        regions[chrom][index] = (rstart, rstop, rlen)

                    mergeCount += 1
                    merged = True
                    break

        if not merged:
            if keepLabel and keepPeak:
                regions[chrom].append((label, start, stop, length, peakPos, peakHeight))
            elif keepLabel:
                regions[chrom].append((label, start, stop, length))
            elif keepPeak:
                regions[chrom].append((start, stop, length, peakPos, peakHeight))
            else:
                regions[chrom].append((start, stop, length))

            count += 1
        if verbose and (count % 100000 == 0):
            print count

    infile.close()

    regionCount = 0
    for chrom in regions:
        regionCount += len(regions[chrom])
        if keepLabel:
            regions[chrom].sort(cmp=lambda x,y:cmp(x[1], y[1]))
        else:
            regions[chrom].sort()

    if verbose:
        print "merged %d times" % mergeCount
        print "returning %d regions" % regionCount

    return regions


def regionsOverlap(start, stop, rstart, rstop):
    return (rstart <= start <= rstop) or (rstart <= stop <= rstop) or (start <= rstart <= stop) or (start <= rstop <= stop)


def regionsAreWithinDistance(start, stop, rstart, rstop, maxDist):
    return (abs(rstart-stop) <= maxDist) or (abs(rstop-start) <= maxDist)


def findPeak(hitList, start, length, readlen=25, doWeight=False, leftPlus=False,
             shift=0, returnShift=False, maxshift=75):
    """ find the peak in a list of reads (hitlist) in a region
    of a given length and absolute start point. returns a
    list of peaks, the number of hits, a triangular-smoothed
    version of hitlist, and the number of reads that are
    forward (plus) sense.
    If doWeight is True, weight the reads accordingly.
    If leftPlus is True, return the number of plus reads left of
    the peak, taken to be the first TopPos position.
    """

    seqArray = array("f", [0.] * length)
    smoothArray = array("f", [0.] * length)
    numHits = 0.
    numPlus = 0.
    regionArray = []
    # TODO: pick best shift for region
    if shift == "auto":
        bestShift = 0
        lowestScore = 20000000000
        for testShift in xrange(maxshift + 1):
            shiftArray = array("f", [0.] * length)
            for read in hitList:
                currentpos = read[0] - start
                if read[1] == "+":
                    currentpos += testShift
                else:
                    currentpos -= testShift

                if (currentpos < 1) or (currentpos >= length):
                    continue

                if doWeight:
                    weight = read[2]
                else:
                    weight = 1.0

                if read[1] == "+":
                    shiftArray[currentpos] += weight
                else:
                    shiftArray[currentpos] -= weight

            currentScore = 0
            for score in shiftArray:
                currentScore += abs(score)

            if currentScore < lowestScore:
                bestShift = testShift
                lowestScore = currentScore

        shift = bestShift

    # once we have the best shift, compute seqArray
    for read in hitList:
        currentpos = read[0] - start
        if read[1] == "+":
            currentpos += shift
        else:
            currentpos -= shift

        if (currentpos <  1 - readlen) or (currentpos >= length):
            continue

        hitIndex = 0
        if doWeight:
            weight = read[2]
        else:
            weight = 1.0

        numHits += weight
        if leftPlus:
            regionArray.append(read)

        while currentpos < 0:
            hitIndex += 1
            currentpos += 1

        while hitIndex < readlen and  currentpos < length:
            seqArray[currentpos] += weight
            hitIndex += 1
            currentpos += 1

        if read[1] == "+":
            numPlus += weight

    # implementing a triangular smooth
    for pos in range(2,length -2):
        smoothArray[pos] = (seqArray[pos -2] + 2 * seqArray[pos - 1] + 3 * seqArray[pos] + 2 * seqArray[pos + 1] + seqArray[pos + 2]) / 9.0

    topNucleotide = 0
    topPos = []
    for currentpos in xrange(length):
        if topNucleotide < smoothArray[currentpos]:
            topNucleotide = smoothArray[currentpos]
            topPos = [currentpos]
        elif topNucleotide  == smoothArray[currentpos]:
            topPos.append(currentpos)

    if leftPlus:
        numLeftPlus = 0
        maxPos = topPos[0]
        for read in regionArray:
            if doWeight:
                weight = read[2]
            else:
                weight = 1.0

            currentPos = read[0] - start
            if currentPos <= maxPos and read[1] == "+":
                numLeftPlus += weight

        if returnShift:
            return (topPos, numHits, smoothArray, numPlus, numLeftPlus, shift)
        else:
            return (topPos, numHits, smoothArray, numPlus, numLeftPlus)
    else:
        if returnShift:
            return (topPos, numHits, smoothArray, numPlus, shift)
        else:
            return (topPos, numHits, smoothArray, numPlus)


def getFeaturesByChromDict(genomeObject, additionalRegionsDict={}, ignorePseudo=False,
                           restrictList=[], regionComplement=False, maxStop=250000000):
    """ return a dictionary of cistematic gene features. Requires
    cistematic, obviously. Can filter-out pseudogenes. Will use
    additional regions dict to supplement gene models, if available.
    Can restrict output to a list of GIDs.
    If regionComplement is set to true, returns the regions *outside* of the
    calculated boundaries, which is useful for retrieving intronic and
    intergenic regions. maxStop is simply used to define the uppermost
    boundary of the complement region.
    """ 
    featuresDict = genomeObject.getallGeneFeatures()
    restrictGID = False
    if len(restrictList) > 0:
        restrictGID = True

    if len(additionalRegionsDict) > 0:
        sortList = []
        for chrom in additionalRegionsDict:
            for (label, start, stop, length) in additionalRegionsDict[chrom]:
                if label not in sortList:
                    sortList.append(label)

                if label not in featuresDict:
                    featuresDict[label] = []
                    sense = "+"
                else:
                    sense = featuresDict[label][0][-1]

                featuresDict[label].append(("custom", chrom, start, stop, sense))

        for gid in sortList:
            featuresDict[gid].sort(cmp=lambda x,y:cmp(x[2], y[2]))

    featuresByChromDict = {}
    for gid in featuresDict:
        if restrictGID and gid not in restrictList:
            continue

        featureList = featuresDict[gid]
        newFeatureList = []
        isPseudo = False
        for (ftype, chrom, start, stop, sense) in featureList:
            if ftype == "PSEUDO":
                isPseudo = True

            if (start, stop, ftype) not in newFeatureList:
                notContained = True
                containedList = []
                for (fstart, fstop, ftype2) in newFeatureList:
                    if start >= fstart and stop <= fstop:
                        notContained = False

                    if start < fstart and stop > fstop:
                        containedList.append((fstart, fstop))

                if len(containedList) > 0:
                    newFList = []
                    notContained = True
                    for (fstart, fstop, ftype2) in newFeatureList:
                        if (fstart, fstop) not in containedList:
                            newFList.append((fstart, fstop, ftype2))
                            if start >= fstart and stop <= fstop:
                                notContained = False

                    newFeatureList = newFList
                if notContained:
                    newFeatureList.append((start, stop, ftype))

        if ignorePseudo and isPseudo:
            continue

        if chrom not in featuresByChromDict:
            featuresByChromDict[chrom] = []

        for (start, stop, ftype) in newFeatureList:
            featuresByChromDict[chrom].append((start, stop, gid, sense, ftype))

    for chrom in featuresByChromDict:
        featuresByChromDict[chrom].sort()

    if regionComplement:
        complementByChromDict = {}
        complementIndex = 0
        for chrom in featuresByChromDict:
            complementByChromDict[chrom] = []
            listLength = len(featuresByChromDict[chrom])
            if listLength > 0:
                currentStart = 0
                for index in range(listLength):
                    currentStop = featuresByChromDict[chrom][index][0]
                    complementIndex += 1
                    if currentStart < currentStop:
                        complementByChromDict[chrom].append((currentStart, currentStop, "nonExon" + str(complementIndex), "F", "nonExon"))

                    currentStart = featuresByChromDict[chrom][index][1]

                currentStop = maxStop
                complementByChromDict[chrom].append((currentStart, currentStop, "nonExon" + str(complementIndex), "F", "nonExon"))

        return (featuresByChromDict, complementByChromDict)
    else:
        return featuresByChromDict


def getLocusByChromDict(genomeObject, upstream = 0, downstream = 0, useCDS=True,
                        additionalRegionsDict={}, ignorePseudo=False, upstreamSpanTSS=False,
                        lengthCDS=0, keepSense=False, adjustToNeighbor=True):
    """ return a dictionary of gene loci. Can be used to retrieve additional
    sequence upstream or downstream of gene, up to the next gene. Requires
    cistematic, obviously.
    Can filter-out pseudogenes and use additional regions outside of existing
    gene models. Use upstreamSpanTSS to overlap half of the upstream region
    over the TSS.
    If lengthCDS > 0 bp, e.g. X, return only the starting X bp from CDS. If
    lengthCDS < 0bp, return only the last X bp from CDS.
    """ 
    locusByChromDict = {}
    if upstream == 0 and downstream == 0 and not useCDS:
        print "getLocusByChromDict: asked for no sequence - returning empty dict"
        return locusByChromDict
    elif upstream > 0 and downstream > 0 and not useCDS:
        print "getLocusByChromDict: asked for only upstream and downstream - returning empty dict"
        return locusByChromDict
    elif lengthCDS != 0 and not useCDS:
        print "getLocusByChromDict: asked for partial CDS but not useCDS - returning empty dict"
        return locusByChromDict
    elif upstreamSpanTSS and lengthCDS != 0:
        print "getLocusByChromDict: asked for TSS spanning and partial CDS - returning empty dict"
        return locusByChromDict
    elif lengthCDS > 0 and downstream > 0:
        print "getLocusByChromDict: asked for discontinuous partial CDS from start and downstream - returning empty dict"
        return locusByChromDict
    elif lengthCDS < 0 and upstream > 0:
        print "getLocusByChromDict: asked for discontinuous partial CDS from stop and upstream - returning empty dict"
        return locusByChromDict

    genome = genomeObject.genome
    featuresDict = genomeObject.getallGeneFeatures()
    if len(additionalRegionsDict) > 0:
        sortList = []
        for chrom in additionalRegionsDict:
            for (label, start, stop, length) in additionalRegionsDict[chrom]:
                if label not in sortList:
                    sortList.append(label)

                if label not in featuresDict:
                    featuresDict[label] = []
                    sense = "+"
                else:
                    sense = featuresDict[label][0][-1]

                featuresDict[label].append(("custom", chrom, start, stop, sense))

        for gid in sortList:
            featuresDict[gid].sort(cmp=lambda x,y:cmp(x[2], y[2]))

    for gid in featuresDict:
        featureList = featuresDict[gid]
        newFeatureList = []
        for (ftype, chrom, start, stop, sense) in featureList:
            newFeatureList.append((start, stop))

        if ignorePseudo and ftype == "PSEUDO":
            continue

        newFeatureList.sort()

        sense = featureList[0][-1]
        gstart = newFeatureList[0][0]
        gstop = newFeatureList[-1][1]
        glen = abs(gstart - gstop)
        if sense == "F":
            if not useCDS and upstream > 0:
                if upstreamSpanTSS:
                    if gstop > (gstart + upstream / 2):
                        gstop = gstart + upstream / 2
                    # otherwise leave gstop as is
                else:
                    gstop = gstart
            elif not useCDS and downstream > 0:
                gstart = gstop

            if upstream > 0:
                if upstreamSpanTSS:
                    distance = upstream / 2
                else:
                    distance = upstream

                if adjustToNeighbor:
                    nextGene = genomeObject.leftGeneDistance((genome, gid), distance * 2)
                    if nextGene < distance * 2:
                        distance = nextGene / 2

                if distance < 1:
                    distance = 1

                gstart -= distance

            if downstream > 0:
                distance = downstream
                if adjustToNeighbor:
                    nextGene = genomeObject.rightGeneDistance((genome, gid), downstream * 2)
                    if nextGene < downstream * 2:
                        distance = nextGene / 2

                if distance < 1:
                    distance = 1

                gstop += distance

            if lengthCDS > 0:
                if lengthCDS < glen:
                    gstop = newFeatureList[0][0] + lengthCDS

            if lengthCDS < 0:
                if abs(lengthCDS) < glen:
                    gstart = newFeatureList[-1][1] + lengthCDS
        else:
            if not useCDS and upstream > 0:
                if upstreamSpanTSS:
                    if gstart < (gstop - upstream / 2):
                        gstart = gstop - upstream / 2
                    # otherwise leave gstart as is
                else:
                    gstart = gstop
            elif not useCDS and downstream > 0:
                    gstop = gstart

            if upstream > 0:
                if upstreamSpanTSS:
                    distance = upstream /2
                else:
                    distance = upstream

                if adjustToNeighbor:
                    nextGene = genomeObject.rightGeneDistance((genome, gid), distance * 2)
                    if nextGene < distance * 2:
                        distance = nextGene / 2

                if distance < 1:
                    distance = 1

                gstop += distance

            if downstream > 0:
                distance = downstream
                if adjustToNeighbor:
                    nextGene = genomeObject.leftGeneDistance((genome, gid), downstream * 2)
                    if nextGene < downstream * 2:
                        distance = nextGene / 2

                if distance < 1:
                    distance = 1

                gstart -= distance

            if lengthCDS > 0:
                if lengthCDS < glen:
                    gstart = newFeatureList[-1][-1] - lengthCDS

            if lengthCDS < 0:
                if abs(lengthCDS) < glen:
                    gstop = newFeatureList[0][0] - lengthCDS

        glen = abs(gstop - gstart)
        if chrom not in locusByChromDict:
            locusByChromDict[chrom] = []

        if keepSense:
            locusByChromDict[chrom].append((gstart, gstop, gid, glen, sense))
        else:
            locusByChromDict[chrom].append((gstart, gstop, gid, glen))

    for chrom in locusByChromDict:
        locusByChromDict[chrom].sort()

    return locusByChromDict


def computeRegionBins(regionsByChromDict, hitDict, bins, readlen, regionList = {},
                      normalizedTag=1., defaultRegionFormat = True, fixedFirstBin=-1,
                      binLength=-1):
    """ returns 2 dictionaries of bin counts and region lengths, given a dictionary of predefined regions,
        a dictionary of reads, a number of bins, the length of reads, and optionally a list of regions
        or a different weight / tag.
    """
    index = 0
    regionsBins = {}
    regionsLen = {}

    if defaultRegionFormat:
        regionIDField = 0
        startField = 1
        stopField = 2
        lengthField = 3
    else:
        startField = 0
        stopField = 1
        regionIDField = 2
        lengthField = 3

    senseField = 4

    print "entering computeRegionBins"
    if len(regionList) > 0:
        for readID in regionList:
            regionsBins[readID] = [0.] * bins
    else:
        for chrom in regionsByChromDict:
            for regionTuple in regionsByChromDict[chrom]:
                regionsBins[regionIDField] = [0.] * bins

    for chrom in hitDict:
        if chrom not in regionsByChromDict:
            continue

        for regionTuple in regionsByChromDict[chrom]:
            regionsLen[regionTuple[regionIDField]] = regionTuple[lengthField]

        print "\n" + chrom
        startRegion = 0
        for (tagStart, sense, weight) in hitDict[chrom]:
            index += 1
            if index % 100000 == 0:
                print "read %d " % index,

            stopPoint = tagStart + readlen
            if startRegion < 0:
                startRegion = 0

            for regionTuple in regionsByChromDict[chrom][startRegion:]:
                start = regionTuple[startField]
                stop = regionTuple[stopField]
                regionID = regionTuple[regionIDField]
                rlen = regionTuple[lengthField]
                try:
                    rsense = regionTuple[senseField]
                except:
                    rsense = "F"

                if tagStart > stop:
                    startRegion += 1
                    continue

                if start > stopPoint:
                    startRegion -= 10
                    break

                if start <= tagStart <= stop:
                    if binLength < 1:
                        regionBinLength = rlen / bins
                    else:
                        regionBinLength = binLength

                    startdist = tagStart - start
                    if rsense == "F":
                        # we are relying on python's integer division quirk
                        binID = startdist / regionBinLength
                        if (fixedFirstBin > 0) and (startdist < fixedFirstBin):
                            binID = 0
                        elif fixedFirstBin > 0:
                            binID = 1

                        if binID >= bins:
                            binID = bins - 1

                        try:
                            regionsBins[regionID][binID] += normalizedTag * weight
                        except:
                            print "%s %s" % (regionID, str(binID))
                    else:
                        rdist = rlen - startdist
                        binID = rdist / regionBinLength
                        if (fixedFirstBin > 0) and (rdist < fixedFirstBin):
                            binID = 0
                        elif fixedFirstBin > 0:
                            binID = 1

                        if binID >= bins:
                            binID = bins - 1

                        try:
                            regionsBins[regionID][binID] += normalizedTag * weight
                        except:
                            print "%s %s" % (regionID, str(binID))

                    stopPoint = stop

    return (regionsBins, regionsLen)


class readDataset:
    """ Class for storing reads from experiments. Assumes that custom scripts
    will translate incoming data into a format that can be inserted into the
    class using the insert* methods. Default class subtype ('DNA') includes
    tables for unique and multireads, whereas 'RNA' subtype also includes a
    splices table.
    """

    def __init__(self, datafile, initialize=False, datasetType='', verbose=False, 
                 cache=False, reportCount=True):
        """ creates an rds datafile if initialize is set to true, otherwise
        will append to existing tables. datasetType can be either 'DNA' or 'RNA'.
        """
        self.dbcon = ""
        self.memcon = ""
        self.dataType = ""
        self.rdsVersion = "1.1"
        self.memBacked = False
        self.memChrom = ""
        self.memCursor = ""
        self.cachedDBFile = ""

        if cache:
            if verbose:
                print "caching ...."

            self.cacheDB(datafile)
            dbfile = self.cachedDBFile
        else:
            dbfile = datafile

        self.dbcon = sqlite.connect(dbfile)
        self.dbcon.row_factory = sqlite.Row
        self.dbcon.execute("PRAGMA temp_store = MEMORY")
        if initialize:
            if datasetType == "":
                self.dataType = "DNA"
            else:
                self.dataType = datasetType

            self.initializeTables(self.dbcon)
        else:
            metadata = self.getMetadata("dataType")
            self.dataType = metadata["dataType"]

        try:
            metadata = self.getMetadata("rdsVersion")
            self.rdsVersion = metadata["rdsVersion"]
        except:
            try:
                self.insertMetadata([("rdsVersion", currentRDSversion)])
            except:
                print "could not add rdsVersion - read-only ?"
                self.rdsVersion = "pre-1.0"

        if verbose:
            if initialize:
                print "INITIALIZED dataset %s" % datafile
            else:
                print "dataset %s" % datafile

            metadata = self.getMetadata()
            print "metadata:"
            pnameList = metadata.keys()
            pnameList.sort()
            for pname in pnameList:
                print "\t" + pname + "\t" + metadata[pname]

            if reportCount:
                ucount = self.getUniqsCount()
                mcount = self.getMultiCount()
                if self.dataType == "DNA" and not initialize:
                    try:
                        print "\n%d unique reads and %d multireads" % (int(ucount), int(mcount))
                    except:
                        print "\n%s unique reads and %s multireads" % (ucount, mcount)
                elif self.dataType == 'RNA' and not initialize:
                    scount = self.getSplicesCount()
                    try:
                        print "\n%d unique reads, %d spliced reads and %d multireads" % (int(ucount), int(scount), int(mcount))
                    except:
                        print "\n%s unique reads, %s spliced reads and %s multireads" % (ucount, scount, mcount)

            print "default cache size is %d pages" % self.getDefaultCacheSize()
            if self.hasIndex():
                print "found index"
            else:
                print "not indexed"


    def __len__(self):
        """ return the number of usable reads in the dataset.
        """
        try:
            total = self.getUniqsCount()
        except:
            total = 0

        try:
            total += self.getMultiCount()
        except:
            pass

        if self.dataType == "RNA":
            try:
                total += self.getSplicesCount()
            except:
                pass

        try:
            total = int(total)
        except:
            total = 0

        return total


    def __del__(self):
        """ cleanup copy in local cache, if present.
        """
        if self.cachedDBFile != "":
            self.uncacheDB()


    def cacheDB(self, filename):
        """ copy geneinfoDB to a local cache.
        """
        self.cachedDBFile = tempfile.mktemp() + ".db"
        shutil.copyfile(filename, self.cachedDBFile)


    def saveCacheDB(self, filename):
        """ copy geneinfoDB to a local cache.
        """
        shutil.copyfile(self.cachedDBFile, filename)


    def uncacheDB(self):
        """ delete geneinfoDB from local cache.
        """
        global cachedDBFile
        if self.cachedDBFile != "":
            try:
                os.remove(self.cachedDBFile)
            except:
                print "could not delete %s" % self.cachedDBFile

            self.cachedDB = ""


    def attachDB(self, filename, asname):
        """ attach another database file to the readDataset.
        """
        stmt = "attach '%s' as %s" % (filename, asname)
        self.execute(stmt)


    def detachDB(self, asname):
        """ detach a database file to the readDataset.
        """
        stmt = "detach %s" % (asname)
        self.execute(stmt)


    def importFromDB(self, asname, table, ascolumns="*", destcolumns="", flagged=""):
        """ import into current RDS the table (with columns destcolumns,
            with default all columns) from the database file asname,
            using the column specification of ascolumns (default all).
        """
        stmt = "insert into %s %s select %s from %s.%s" % (table, destcolumns, ascolumns, asname, table)
        if flagged != "":
            stmt += " where flag = '%s' " % flagged

        self.execute(stmt, forceCommit=True)


    def getTables(self, asname=""):
        """ get a list of table names in a particular database file.
        """
        resultList = []

        if self.memBacked:
            sql = self.memcon.cursor()
        else:
            sql = self.dbcon.cursor()

        if asname != "":
            asname += "."

        stmt = "select name from %ssqlite_master where type='table'" % asname
        sql.execute(stmt)
        results = sql.fetchall()

        for row in results:
            resultList.append(row["name"])

        return resultList


    def hasIndex(self):
        """ check whether the RDS file has at least one index.
        """
        stmt = "select count(*) from sqlite_master where type='index'"
        count = int(self.execute(stmt, returnResults=True)[0][0])
        if count > 0:
            return True

        return False


    def initializeTables(self, acon, cache=100000):
        """ creates table schema in database connection acon, which is
        typically a database file or an in-memory database.
        """
        acon.execute("PRAGMA DEFAULT_CACHE_SIZE = %d" % cache)
        acon.execute("create table metadata (name varchar, value varchar)")
        acon.execute("insert into metadata values('dataType','%s')" % self.dataType)
        acon.execute("create table uniqs (ID INTEGER PRIMARY KEY, readID varchar, chrom varchar, start int, stop int, sense varchar, weight real, flag varchar, mismatch varchar)")
        acon.execute("create table multi (ID INTEGER PRIMARY KEY, readID varchar, chrom varchar, start int, stop int, sense varchar, weight real, flag varchar, mismatch varchar)")
        if self.dataType == "RNA":
            acon.execute("create table splices (ID INTEGER PRIMARY KEY, readID varchar, chrom varchar, startL int, stopL int, startR int, stopR int, sense varchar, weight real, flag varchar, mismatch varchar)")

        acon.commit()


    def getFileCursor(self):
        """ returns a cursor to file database for low-level (SQL)
        access to the data.
        """
        return self.dbcon.cursor()


    def getMemCursor(self):
        """ returns a cursor to memory database for low-level (SQL)
        access to the data.
        """
        return self.memcon.cursor()


    def getMetadata(self, valueName=""):
        """ returns a dictionary of metadata.
        """
        whereClause = ""
        resultsDict = {}

        if valueName != "":
            whereClause = " where name = '%s' " % valueName

        if self.memBacked:
            sql = self.memcon.cursor()
        else:
            sql = self.dbcon.cursor()

        sql.execute("select name, value from metadata" + whereClause)
        results = sql.fetchall()

        for row in results:
            pname = row["name"]
            pvalue = row["value"]
            if pname not in resultsDict:
                resultsDict[pname] = pvalue
            else:
                trying = True
                index = 2
                while trying:
                    newName = pname + ":" + str(index)
                    if newName not in resultsDict:
                        resultsDict[newName] = pvalue
                        trying = False

                    index += 1

        return resultsDict


    def getReadSize(self):
        """ returns readsize if defined in metadata.
        """
        metadata = self.getMetadata()
        if "readsize" not in metadata:
            print "no readsize parameter defined - returning 0"
            return 0
        else:
            mysize = metadata["readsize"]
            if "import" in mysize:
                mysize = mysize.split()[0]

            return int(mysize)


    def getDefaultCacheSize(self):
        """ returns the default cache size.
        """
        return int(self.execute("PRAGMA DEFAULT_CACHE_SIZE", returnResults=True)[0][0])


    def getChromosomes(self, table="uniqs", fullChrom=True):
        """ returns a list of distinct chromosomes in table.
        """
        statement = "select distinct chrom from %s" % table
        if self.memBacked:
            sql = self.memcon.cursor()
        else:
            sql = self.dbcon.cursor()

        sql.execute(statement)
        results = []
        for row in sql:
            if fullChrom:
                if row["chrom"] not in results:
                    results.append(row["chrom"])
            else:
                if  len(row["chrom"][3:].strip()) < 1:
                    continue

                if row["chrom"][3:] not in results:
                    results.append(row["chrom"][3:])

        results.sort()

        return results


    def getMaxCoordinate(self, chrom, verbose=False, doUniqs=True,
                         doMulti=False, doSplices=False):
        """ returns the maximum coordinate for reads on a given chromosome.
        """
        maxCoord = 0
        if self.memBacked:
            sql = self.memcon.cursor()
        else:
            sql = self.dbcon.cursor()

        if doUniqs:
            try:
                sql.execute("select max(start) from uniqs where chrom = '%s'" % chrom)
                maxCoord = int(sql.fetchall()[0][0])
            except:
                print "couldn't retrieve coordMax for chromosome %s" % chrom

        if doSplices:
            sql.execute("select max(startR) from splices where chrom = '%s'" % chrom)
            try:
                spliceMax = int(sql.fetchall()[0][0])
                if spliceMax > maxCoord:
                    maxCoord = spliceMax
            except:
                pass

        if doMulti:
            sql.execute("select max(start) from multi where chrom = '%s'" % chrom)
            try:
                multiMax = int(sql.fetchall()[0][0])
                if multiMax > maxCoord:
                    maxCoord = multiMax
            except:
                pass

        if verbose:
            print "%s maxCoord: %d" % (chrom, maxCoord)

        return maxCoord


    def getReadsDict(self, verbose=False, bothEnds=False, noSense=False, fullChrom=False, chrom="",
                     flag="", withWeight=False, withFlag=False, withMismatch=False, withID=False,
                     withChrom=False, withPairID=False, doUniqs=True, doMulti=False, findallOptimize=False,
                     readIDDict=False, readLike="", start=-1, stop=-1, limit=-1, hasMismatch=False,
                     flagLike=False, strand="", entryDict=False):
        """ returns a dictionary of reads in a variety of formats
        and which can be restricted by chromosome or custom-flag.
        Returns unique reads by default, but can return multireads
        with doMulti set to True.
        """
        whereClause = []
        resultsDict = {}

        if chrom != "" and chrom != self.memChrom:
            whereClause.append("chrom = '%s'" % chrom)

        if flag != "":
            if flagLike:
                flagLikeClause = string.join(['flag LIKE "%', flag, '%"'], "")
                whereClause.append(flagLikeClause)
            else:
                whereClause.append("flag = '%s'" % flag)

        if start > -1:
            whereClause.append("start > %d" % start)

        if stop > -1:
            whereClause.append("stop < %d" % stop)

        if len(readLike) > 0:
            readIDClause = string.join(["readID LIKE  '", readLike, "%'"], "")
            whereClause.append(readIDClause)

        if hasMismatch:
            whereClause.append("mismatch != ''")

        if strand in ["+", "-"]:
            whereClause.append("sense = '%s'" % strand)

        if len(whereClause) > 0:
            whereStatement = string.join(whereClause, " and ")
            whereQuery = "where %s" % whereStatement
        else:
            whereQuery = ""

        groupBy = []
        if findallOptimize:
            selectClause = ["select start, sense, sum(weight)"]
            groupBy = ["GROUP BY start, sense"]
        else:
            selectClause = ["select ID, chrom, start, readID"]
            if bothEnds:
                selectClause.append("stop")

            if not noSense:
                selectClause.append("sense")

            if withWeight:
                selectClause.append("weight")

            if withFlag:
                selectClause.append("flag")

            if withMismatch:
                selectClause.append("mismatch")

        if limit > 0:
            groupBy.append("LIMIT %d" % limit)

        selectQuery = string.join(selectClause, ",")
        groupQuery = string.join(groupBy)
        if doUniqs:
            stmt = [selectQuery, "from uniqs", whereQuery, groupQuery]
            if doMulti:
                stmt.append("UNION ALL")
                stmt.append(selectQuery)
                stmt.append("from multi")
                stmt.append(whereQuery)
                stmt.append(groupQuery)
        else:
            stmt = [selectQuery, "from multi", whereQuery]

        if findallOptimize:
            if self.memBacked:
                self.memcon.row_factory = None
                sql = self.memcon.cursor()
            else:
                self.dbcon.row_factory = None
                sql = self.dbcon.cursor()

            stmt.append("order by start")
        elif readIDDict:
            if self.memBacked:
                sql = self.memcon.cursor()
            else:
                sql = self.dbcon.cursor()

            stmt.append("order by readID, start")
        else:
            if self.memBacked:
                sql = self.memcon.cursor()
            else:
                sql = self.dbcon.cursor()

            stmt.append("order by chrom, start")

        sqlQuery = string.join(stmt)
        sql.execute(sqlQuery)

        if findallOptimize:
            resultsDict[chrom] = [[int(row[0]), row[1], float(row[2])] for row in sql]
            if self.memBacked:
                self.memcon.row_factory = sqlite.Row
            else:
                self.dbcon.row_factory = sqlite.Row
        else:
            currentChrom = ""
            currentReadID = ""
            pairID = 0
            for row in sql:
                readID = row["readID"]
                if fullChrom:
                    chrom = row["chrom"]
                else:
                    chrom = row["chrom"][3:]

                if not readIDDict and chrom != currentChrom:
                    resultsDict[chrom] = []
                    currentChrom = chrom
                    dictKey = chrom
                elif readIDDict:
                    theReadID = readID
                    if "::" in readID:
                        (theReadID, multiplicity) = readID.split("::")

                    if "/" in theReadID and withPairID:
                        (theReadID, pairID) = readID.split("/")

                    if theReadID != currentReadID:
                        resultsDict[theReadID] = []
                        currentReadID = theReadID
                        dictKey = theReadID

                if entryDict:
                    newrow = {"start": int(row["start"])}
                    if bothEnds:
                        newrow["stop"] = int(row["stop"])

                    if not noSense:
                        newrow["sense"] = row["sense"]

                    if withWeight:
                        newrow["weight"] = float(row["weight"])

                    if withFlag:
                        newrow["flag"] = row["flag"]

                    if withMismatch:
                        newrow["mismatch"] = row["mismatch"]

                    if withID:
                        newrow["readID"] = readID

                    if withChrom:
                        newrow["chrom"] = chrom

                    if withPairID:
                        newrow["pairID"] = pairID
                else:
                    newrow = [int(row["start"])]
                    if bothEnds:
                        newrow.append(int(row["stop"]))

                    if not noSense:
                        newrow.append(row["sense"])

                    if withWeight:
                        newrow.append(float(row["weight"]))

                    if withFlag:
                        newrow.append(row["flag"])

                    if withMismatch:
                        newrow.append(row["mismatch"])

                    if withID:
                        newrow.append(readID)

                    if withChrom:
                        newrow.append(chrom)

                    if withPairID:
                        newrow.append(pairID)

                resultsDict[dictKey].append(newrow)

        return resultsDict


    def getSplicesDict(self, verbose=False, noSense=False, fullChrom=False, chrom="",
                       flag="", withWeight=False, withFlag=False, withMismatch=False,
                       withID=False, withChrom=False, withPairID=False, readIDDict=False,
                       splitRead=False, hasMismatch=False, flagLike=False, start=-1,
                       stop=-1, strand="", entryDict=False):
        """ returns a dictionary of spliced reads in a variety of
        formats and which can be restricted by chromosome or custom-flag.
        Returns unique spliced reads for now.
        """
        whereClause = []
        resultsDict = {}

        if chrom != "" and chrom != self.memChrom:
            whereClause = ["chrom = '%s'" % chrom]

        if flag != "":
            if flagLike:
                flagLikeClause = string.join(['flag LIKE "%', flag, '%"'], "")
                whereClause.append(flagLikeClause)
            else:
                whereClause.append("flag = '%s'" % flag)

        if hasMismatch:
            whereClause.append("mismatch != ''")

        if strand != "":
            whereClause.append("sense = '%s'" % strand)

        if start > -1:
            whereClause.append("startL > %d" % start)

        if stop > -1:
            whereClause.append("stopR < %d" % stop)

        if len(whereClause) > 0:
            whereStatement = string.join(whereClause, " and ")
            whereQuery = "where %s" % whereStatement
        else:
            whereQuery = ""

        selectClause = ["select ID, chrom, startL, stopL, startR, stopR, readID"]
        if not noSense:
            selectClause.append("sense")

        if withWeight:
            selectClause.append("weight")

        if withFlag:
            selectClause.append("flag")

        if withMismatch:
            selectClause.append("mismatch")

        selectQuery = string.join(selectClause, " ,")
        if self.memBacked:
            sql = self.memcon.cursor()
        else:
            sql = self.dbcon.cursor()

        if chrom == "" and not readIDDict:
            stmt = "select distinct chrom from splices %s" % whereQuery
            sql.execute(stmt)
            for row in sql:
                if fullChrom:
                    chrom = row["chrom"]
                else:
                    chrom = row["chrom"][3:]

                resultsDict[chrom] = []
        elif chrom != "" and not readIDDict:
            resultsDict[chrom] = []

        stmt = "%s from splices %s order by chrom, startL" % (selectQuery, whereQuery)
        sql.execute(stmt)
        currentReadID = ""
        for row in sql:
            pairID = 0
            readID = row["readID"]
            if fullChrom:
                chrom = row["chrom"]
            else:
                chrom = row["chrom"][3:]

            if readIDDict:
                if "/" in readID:
                    (theReadID, pairID) = readID.split("/")
                else:
                    theReadID = readID

                if theReadID != currentReadID:
                    resultsDict[theReadID] = []
                    currentReadID = theReadID
                    dictKey = theReadID
            else:
                dictKey = chrom

            if entryDict:
                newrow = {"startL": int(row["startL"])}
                newrow["stopL"] = int(row["stopL"])
                newrow["startR"] = int(row["startR"])
                newrow["stopR"] = int(row["stopR"])
                if not noSense:
                    newrow["sense"] = row["sense"]

                if withWeight:
                    newrow["weight"] = float(row["weight"])

                if withFlag:
                    newrow["flag"] = row["flag"]

                if withMismatch:
                    newrow["mismatch"] = row["mismatch"]

                if withID:
                    newrow["readID"] = readID

                if withChrom:
                    newrow["chrom"] = chrom

                if withPairID:
                    newrow["pairID"] = pairID

                if splitRead:
                    leftDict = newrow
                    del leftDict["startR"]
                    del leftDict["stopR"]
                    rightDict = newrow
                    del rightDict["start"]
                    del rightDict["stopL"]
                    resultsDict[dictKey].append(leftDict)
                    resultsDict[dictKey].append(rightDict)
                else:
                    resultsDict[dictKey].append(newrow)
            else:
                newrow = [int(row["startL"])]
                newrow.append(int(row["stopL"]))
                newrow.append(int(row["startR"]))
                newrow.append(int(row["stopR"]))
                if not noSense:
                    newrow.append(row["sense"])

                if withWeight:
                    newrow.append(float(row["weight"]))

                if withFlag:
                    newrow.append(row["flag"])

                if withMismatch:
                    newrow.append(row["mismatch"])

                if withID:
                    newrow.append(readID)

                if withChrom:
                    newrow.append(chrom)

                if withPairID:
                    newrow.append(pairID)

                if splitRead:
                    resultsDict[dictKey].append(newrow[:2] + newrow[4:])
                    resultsDict[dictKey].append(newrow[2:])
                else:
                    resultsDict[dictKey].append(newrow)

        return resultsDict


    def getCounts(self, chrom="", rmin="", rmax="", uniqs=True, multi=False,
                  splices=False, reportCombined=True, sense="both"):
        """ return read counts for a given region.
        """
        ucount = 0
        mcount = 0
        scount = 0
        restrict = ""
        if sense in ["+", "-"]:
            restrict = " sense ='%s' " % sense

        if uniqs:
            try:
                ucount = float(self.getUniqsCount(chrom, rmin, rmax, restrict))
            except:
                ucount = 0

        if multi:
            try:
                mcount = float(self.getMultiCount(chrom, rmin, rmax, restrict))
            except:
                mcount = 0

        if splices:
            try:
                scount = float(self.getSplicesCount(chrom, rmin, rmax, restrict))
            except:
                scount = 0

        if reportCombined:
            total = ucount + mcount + scount
            return total
        else:
            return (ucount, mcount, scount)


    def getTotalCounts(self, chrom="", rmin="", rmax=""):
        """ return read counts for a given region.
        """
        return self.getCounts(chrom, rmin, rmax, uniqs=True, multi=True, splices=True, reportCombined=True, sense="both")


    def getTableEntryCount(self, table, chrom="", rmin="", rmax="", restrict="", distinct=False, startField="start"):
        """ returns the number of row in the uniqs table.
        """
        whereClause = []
        count = 0

        if chrom !=""  and chrom != self.memChrom:
            whereClause = ["chrom='%s'" % chrom]

        if rmin != "":
            whereClause.append("%s >= %s" % (startField, str(rmin)))

        if rmax != "":
            whereClause.append("%s <= %s" % (startField, str(rmax)))

        if restrict != "":
            whereClause.append(restrict)

        if len(whereClause) > 0:
            whereStatement = string.join(whereClause, " and ")
            whereQuery = "where %s" % whereStatement
        else:
            whereQuery = ""

        if self.memBacked:
            sql = self.memcon.cursor()
        else:
            sql = self.dbcon.cursor()

        if distinct:
            sql.execute("select count(distinct chrom+start+sense) from %s%s" % (table, whereQuery))
        else:
            sql.execute("select sum(weight) from %s%s" % (table, whereQuery))

        result = sql.fetchone()

        try:
            count = int(result[0])
        except:
            count = 0

        return count


    def getSplicesCount(self, chrom="", rmin="", rmax="", restrict="", distinct=False):
        """ returns the number of row in the splices table.
        """
        return self.getTableEntryCount("splices", chrom, rmin, rmax, restrict, distinct, startField="startL")


    def getUniqsCount(self, chrom="", rmin="", rmax="", restrict="", distinct=False):
        """ returns the number of distinct readIDs in the uniqs table.
        """
        return self.getTableEntryCount("uniqs", chrom, rmin, rmax, restrict, distinct)


    def getMultiCount(self, chrom="", rmin="", rmax="", restrict="", distinct=False):
        """ returns the total weight of readIDs in the multi table.
        """
        return self.getTableEntryCount("multi", chrom, rmin, rmax, restrict, distinct)


    def getReadIDs(self, uniqs=True, multi=False, splices=False, paired=False, limit=-1):
        """ get readID's.
        """
        stmt = []
        limitPart = ""
        if limit > 0:
            limitPart = "LIMIT %d" % limit

        if uniqs:
            stmt.append("select readID from uniqs")

        if multi:
            stmt.append("select readID from multi")

        if splices:
            stmt.append("select readID from splices")

        if len(stmt) > 0:
            selectPart = string.join(stmt, " union ")
        else:
            selectPart = ""

        sqlQuery = "%s group by readID %s" (selectPart, limitPart)
        if self.memBacked:
            sql = self.memcon.cursor()
        else:
            sql = self.dbcon.cursor()

        sql.execute(sqlQuery)
        result = sql.fetchall()

        if paired:
            return [x.split("/")[0][0] for x in result]
        else:
            return [x[0] for x in result]


    def getMismatches(self, mischrom = None, verbose=False, useSplices=True):
        """ returns the uniq and spliced mismatches in a dictionary.
        """
        revcomp = {"A": "T",
                   "T": "A",
                   "G": "C",
                   "C": "G",
                   "N": "N"
        }

        readlen = self.getReadSize()
        if mischrom:
            hitChromList = [mischrom]
        else:
            hitChromList = self.getChromosomes()
            hitChromList.sort()

        snpDict = {}
        for achrom in hitChromList:
            if verbose:
                print "getting mismatches from chromosome %s" % (achrom)

            snpDict[achrom] = []
            hitDict = self.getReadsDict(fullChrom=True, chrom=achrom, withMismatch=True, findallOptimize=False, hasMismatch=True)
            if useSplices and self.dataType == "RNA":
                spliceDict = self.getSplicesDict(fullChrom=True, chrom=achrom, withMismatch=True, readIDDict=True, hasMismatch=True)
                spliceIDList = spliceDict.keys()
                for k in spliceIDList:
                    (startpos, lefthalf, rightstart, endspos, sense, mismatches) = spliceDict[k][0]
                    spMismatchList = mismatches.split(",")
                    for mismatch in spMismatchList:
                        if "N" in mismatch:
                            continue

                        change_len = len(mismatch)
                        if sense == "+":
                            change_from = mismatch[0]
                            change_base = mismatch[change_len-1]
                            change_pos = int(mismatch[1:change_len-1])
                        elif sense == "-":
                            change_from = revcomp[mismatch[0]]
                            change_base = revcomp[mismatch[change_len-1]]
                            change_pos = readlen - int(mismatch[1:change_len-1]) + 1

                        firsthalf = int(lefthalf)-int(startpos)+1
                        secondhalf = 0
                        if int(change_pos) <= int(firsthalf):
                            change_at = startpos + change_pos - 1
                        else:
                            secondhalf = change_pos - firsthalf
                            change_at = rightstart + secondhalf

                        snpDict[achrom].append([startpos, change_at, change_base, change_from])

            if achrom not in hitDict:
                continue

            for (start, sense, mismatches) in hitDict[achrom]:
                mismatchList = mismatches.split(",")
                for mismatch in mismatchList:
                    if "N" in mismatch:
                        continue

                    change_len = len(mismatch)
                    if sense == "+":
                        change_from = mismatch[0]
                        change_base = mismatch[change_len-1]
                        change_pos = int(mismatch[1:change_len-1])
                    elif sense == "-":
                        change_from = revcomp[mismatch[0]]
                        change_base = revcomp[mismatch[change_len-1]]
                        change_pos = readlen - int(mismatch[1:change_len-1]) + 1

                    change_at = start + change_pos - 1
                    snpDict[achrom].append([start, change_at, change_base, change_from])

        return snpDict


    def getChromProfile(self, chromosome, cstart=-1, cstop=-1, useMulti=True,
                        useSplices=False, normalizationFactor = 1.0, trackStrand=False,
                        keepStrand="both", shiftValue=0):
        """return a profile of the chromosome as an array of per-base read coverage....
            keepStrand = 'both', 'plusOnly', or 'minusOnly'.
            Will also shift position of unique and multireads (but not splices) if shift is a natural number
        """
        metadata = self.getMetadata()
        readlen = int(metadata["readsize"])
        dataType = metadata["dataType"]
        scale = 1. / normalizationFactor
        shift = {}
        shift["+"] = int(shiftValue)
        shift["-"] = -1 * int(shiftValue)

        if cstop > 0:
            lastNT = self.getMaxCoordinate(chromosome, doMulti=useMulti, doSplices=useSplices) + readlen
        else:
            lastNT = cstop - cstart + readlen + shift["+"]

        chromModel = array("f", [0.] * lastNT)
        hitDict = self.getReadsDict(fullChrom=True, chrom=chromosome, withWeight=True, doMulti=useMulti, start=cstart, stop=cstop, findallOptimize=True)
        if cstart < 0:
            cstart = 0

        for (hstart, sense, weight) in hitDict[chromosome]:
            hstart = hstart - cstart + shift[sense]
            for currentpos in range(hstart,hstart+readlen):
                try:
                    if not trackStrand or (sense == "+" and keepStrand != "minusOnly"):
                        chromModel[currentpos] += scale * weight
                    elif sense == '-' and keepStrand != "plusOnly":
                        chromModel[currentpos] -= scale * weight
                except:
                    continue

        del hitDict
        if useSplices and dataType == "RNA":
            if cstop > 0:
                spliceDict = self.getSplicesDict(fullChrom=True, chrom=chromosome, withID=True, start=cstart, stop=cstop)
            else:
                spliceDict = self.getSplicesDict(fullChrom=True, chrom=chromosome, withID=True)
   
            if chromosome in spliceDict:
                for (Lstart, Lstop, Rstart, Rstop, rsense, readName) in spliceDict[chromosome]:
                    if (Rstop - cstart) < lastNT:
                        for index in range(abs(Lstop - Lstart)):
                            currentpos = Lstart - cstart + index
                            # we only track unique splices
                            if not trackStrand or (rsense == "+" and keepStrand != "minusOnly"):
                                chromModel[currentpos] += scale
                            elif rsense == "-" and keepStrand != "plusOnly":
                                chromModel[currentpos] -= scale

                        for index in range(abs(Rstop - Rstart)):
                            currentpos = Rstart - cstart + index
                            # we only track unique splices
                            if not trackStrand or (rsense == "+" and keepStrand != "minusOnly"):
                                chromModel[currentpos] += scale
                            elif rsense == "-" and keepStrand != "plusOnly":
                                chromModel[currentpos] -= scale

            del spliceDict

        return chromModel


    def insertMetadata(self, valuesList):
        """ inserts a list of (pname, pvalue) into the metadata
        table.
        """
        self.dbcon.executemany("insert into metadata(name, value) values (?,?)", valuesList)
        self.dbcon.commit()


    def updateMetadata(self, pname, newValue, originalValue=""):
        """ update a metadata field given the original value and the new value.
        """
        stmt = "update metadata set value='%s' where name='%s'" % (str(newValue), pname)
        if originalValue != "":
            stmt += " and value='%s' " % str(originalValue)

        self.dbcon.execute(stmt)
        self.dbcon.commit()


    def insertUniqs(self, valuesList):
        """ inserts a list of (readID, chrom, start, stop, sense, weight, flag, mismatch)
        into the uniqs table.
        """
        self.dbcon.executemany("insert into uniqs(ID, readID, chrom, start, stop, sense, weight, flag, mismatch) values (NULL,?,?,?,?,?,?,?,?)", valuesList)
        self.dbcon.commit()


    def insertMulti(self, valuesList):
        """ inserts a list of (readID, chrom, start, stop, sense, weight, flag, mismatch)
        into the multi table.
        """
        self.dbcon.executemany("insert into multi(ID, readID, chrom, start, stop, sense, weight, flag, mismatch) values (NULL,?,?,?,?,?,?,?,?)", valuesList)
        self.dbcon.commit()


    def insertSplices(self, valuesList):
        """ inserts a list of (readID, chrom, startL, stopL, startR, stopR, sense, weight, flag, mismatch)
        into the splices table.
        """
        self.dbcon.executemany("insert into splices(ID, readID, chrom, startL, stopL, startR, stopR, sense, weight, flag, mismatch) values (NULL,?,?,?,?,?,?,?,?,?,?)", valuesList)
        self.dbcon.commit()


    def flagReads(self, regionsList, uniqs=True, multi=False, splices=False, sense="both"):
        """ update reads on file database in a list region of regions for a chromosome to have a new flag.
            regionsList must have 4 fields per region of the form (flag, chrom, start, stop) or, with
            sense set to '+' or '-', 5 fields per region of the form (flag, chrom, start, stop, sense).
        """
        restrict = ""
        if sense != "both":
            restrict = " and sense = ? "

        if uniqs:
            self.dbcon.executemany("UPDATE uniqs SET flag = ? where chrom = ? and start >= ? and start < ? " + restrict, regionsList)

        if multi:
            self.dbcon.executemany("UPDATE multi SET flag = ? where chrom = ? and start >= ? and start < ? " + restrict, regionsList)

        if self.dataType == "RNA" and splices:
            self.dbcon.executemany("UPDATE splices SET flag = flag || ' L:' || ? where chrom = ? and startL >= ? and startL < ? " + restrict, regionsList)
            self.dbcon.executemany("UPDATE splices SET flag = flag || ' R:' || ? where chrom = ? and startR >= ? and startR < ? " + restrict, regionsList)

        self.dbcon.commit()


    def setFlags(self, flag, uniqs=True, multi=True, splices=True):
        """ set the flag fields in the entire dataset to clear. Useful for rerunning an analysis from scratch.
        """
        if uniqs:
            self.dbcon.execute("UPDATE uniqs SET flag = '%s'" % flag)

        if multi:
            self.dbcon.execute("UPDATE multi SET flag = '%s'" % flag)

        if self.dataType == 'RNA' and splices:
            self.dbcon.execute("UPDATE splices SET flag = '%s'" % flag)

        self.dbcon.commit()


    def resetFlags(self, uniqs=True, multi=True, splices=True):
        """ reset the flag fields in the entire dataset to clear. Useful for rerunning an analysis from scratch.
        """
        if uniqs:
            self.dbcon.execute("UPDATE uniqs SET flag = ''")

        if multi:
            self.dbcon.execute("UPDATE multi SET flag = ''")

        if self.dataType == "RNA" and splices:
            self.dbcon.execute("UPDATE splices SET flag = ''")

        self.dbcon.commit()


    def reweighMultireads(self, readList):
        self.dbcon.executemany("UPDATE multi SET weight = ? where chrom = ? and start = ? and readID = ? ", readList)


    def setSynchronousPragma(self, value="ON"):
        try:
            self.dbcon.execute("PRAGMA SYNCHRONOUS = %s" % value)
        except:
            print "warning: couldn't set PRAGMA SYNCHRONOUS = %s" % value


    def setDBcache(self, cache, default=False):
        self.dbcon.execute("PRAGMA CACHE_SIZE = %d" % cache)
        if default:
            self.dbcon.execute('PRAGMA DEFAULT_CACHE_SIZE = %d' % cache)


    def execute(self, statement, returnResults=False, forceCommit=False):
        if self.memBacked:
            sql = self.memcon.cursor()
        else:
            sql = self.dbcon.cursor()

        sql.execute(statement)
        if returnResults:
            result = sql.fetchall()
            return result

        if forceCommit:
            if self.memBacked:
                self.memcon.commit()
            else:
                self.dbcon.commit()


    def buildIndex(self, cache=100000):
        """ Builds the file indeces for the main tables.
            Cache is the number of 1.5 kb pages to keep in memory.
            100000 pages translates into 150MB of RAM, which is our default.
        """
        if cache > self.getDefaultCacheSize():
            self.setDBcache(cache)
        self.setSynchronousPragma("OFF")
        self.dbcon.execute("CREATE INDEX uPosIndex on uniqs(chrom, start)")
        print "built uPosIndex"
        self.dbcon.execute("CREATE INDEX uChromIndex on uniqs(chrom)")
        print "built uChromIndex"
        self.dbcon.execute("CREATE INDEX mPosIndex on multi(chrom, start)")
        print "built mPosIndex"
        self.dbcon.execute("CREATE INDEX mChromIndex on multi(chrom)")
        print "built mChromIndex"

        if self.dataType == "RNA":
            self.dbcon.execute("CREATE INDEX sPosIndex on splices(chrom, startL)")
            print "built sPosIndex"
            self.dbcon.execute("CREATE INDEX sPosIndex2 on splices(chrom, startR)")
            print "built sPosIndex2"
            self.dbcon.execute("CREATE INDEX sChromIndex on splices(chrom)")
            print "built sChromIndex"

        self.dbcon.commit()
        self.setSynchronousPragma("ON")


    def dropIndex(self):
        """ drops the file indices for the main tables.
        """
        try:
            self.setSynchronousPragma("OFF")
            self.dbcon.execute("DROP INDEX uPosIndex")
            self.dbcon.execute("DROP INDEX uChromIndex")
            self.dbcon.execute("DROP INDEX mPosIndex")
            self.dbcon.execute("DROP INDEX mChromIndex")

            if self.dataType == "RNA":
                self.dbcon.execute("DROP INDEX sPosIndex")
                try:
                    self.dbcon.execute("DROP INDEX sPosIndex2")
                except:
                    pass

                self.dbcon.execute("DROP INDEX sChromIndex")

            self.dbcon.commit()
        except:
            print "problem dropping index"

        self.setSynchronousPragma("ON")


    def memSync(self, chrom="", index=False):
        """ makes a copy of the dataset into memory for faster access.
        Can be restricted to a "full" chromosome. Can also build the
        memory indices.
        """
        self.memcon = ""
        self.memcon = sqlite.connect(":memory:")
        self.initializeTables(self.memcon)
        cursor = self.dbcon.cursor()
        whereclause = ""
        if chrom != "":
            print "memSync %s" % chrom
            whereclause = " where chrom = '%s' " % chrom
            self.memChrom = chrom
        else:
            self.memChrom = ""

        self.memcon.execute("PRAGMA temp_store = MEMORY")
        self.memcon.execute("PRAGMA CACHE_SIZE = 1000000")
        # copy metadata to memory
        self.memcon.execute("delete from metadata")
        results = cursor.execute("select name, value from metadata")
        results2 = []
        for row in results:
            results2.append((row["name"], row["value"]))

        self.memcon.executemany("insert into metadata(name, value) values (?,?)", results2)
        # copy uniqs to memory
        results = cursor.execute("select chrom, start, stop, sense, weight, flag, mismatch, readID from uniqs" + whereclause)
        results2 = []
        for row in results:
            results2.append((row["readID"], row["chrom"], int(row["start"]), int(row["stop"]), row["sense"], row["weight"], row["flag"], row["mismatch"]))

        self.memcon.executemany("insert into uniqs(ID, readID, chrom, start, stop, sense, weight, flag, mismatch) values (NULL,?,?,?,?,?,?,?,?)", results2)
        # copy multi to memory
        results = cursor.execute("select chrom, start, stop, sense, weight, flag, mismatch, readID from multi" + whereclause)
        results2 = []
        for row in results:
            results2.append((row["readID"], row["chrom"], int(row["start"]), int(row["stop"]), row["sense"], row["weight"], row["flag"], row["mismatch"]))

        self.memcon.executemany("insert into multi(ID, readID, chrom, start, stop, sense, weight, flag, mismatch) values (NULL,?,?,?,?,?,?,?,?)", results2)
        # copy splices to memory
        if self.dataType == "RNA":
            results = cursor.execute("select chrom, startL, stopL, startR, stopR, sense, weight, flag, mismatch, readID from splices" + whereclause)
            results2 = []
            for row in results:
                results2.append((row["readID"], row["chrom"], int(row["startL"]), int(row["stopL"]), int(row["startR"]), int(row["stopR"]), row["sense"], row["weight"], row["flag"], row["mismatch"]))

            self.memcon.executemany("insert into splices(ID, readID, chrom, startL, stopL, startR, stopR, weight, sense, flag, mismatch) values (NULL,?,?,?,?,?,?,?,?,?,?)", results2)
        if index:
            if chrom != "":
                self.memcon.execute("CREATE INDEX uPosIndex on uniqs(start)")
                self.memcon.execute("CREATE INDEX mPosIndex on multi(start)")
                if self.dataType == "RNA":
                    self.memcon.execute("CREATE INDEX sPosLIndex on splices(startL)")
                    self.memcon.execute("CREATE INDEX sPosRIndex on splices(startR)")
            else:
                self.memcon.execute("CREATE INDEX uPosIndex on uniqs(chrom, start)")
                self.memcon.execute("CREATE INDEX mPosIndex on multi(chrom, start)")
                if self.dataType == "RNA":
                    self.memcon.execute("CREATE INDEX sPosLIndex on splices(chrom, startL)")
                    self.memcon.execute("CREATE INDEX sPosRIndex on splices(chrom, startR)")

        self.memBacked = True
        self.memcon.row_factory = sqlite.Row
        self.memcon.commit()
