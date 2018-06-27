import sys, array
from numpy import *

versionString = '%s: version 0.95' % sys.argv[0]
print versionString

def compNT(nt):
    """ returns the complementary basepair to base nt
    """
    compDict = { 'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'S': 'S', 'W': 'W', 'R': 'Y', 'Y': 'R',
                 'M': 'K', 'K': 'M', 'H': 'D', 'D': 'H', 'B': 'V', 'V': 'B', 'N': 'N', 'a': 't',
                 't': 'a', 'g': 'c', 'c': 'g', 'n': 'n' , 'z': 'z'}
    return compDict.get(nt, 'N')

def complement(sequence, length=-1):
    """ returns the complement of the sequence.
    """
    newSeq = ""
    
    seqLength = len(sequence)
    
    if length == seqLength or length < 0:
        seqList = list(sequence)
        seqList.reverse()
        return "".join(map(compNT, seqList))
    
    for index in range(seqLength - 1,seqLength - length - 1, -1):
        try:
            newSeq += compNT(sequence[index])
        except:
            newSeq += 'N'
    #print "length of seq is %d and newseq is %d\n" % (len(sequence), len(newSeq))
    return newSeq

if len(sys.argv) < 5:
    print "python %s incontigfile distalPairs outpathfile outcontigfile [-prefix string] [-overlap bp]\n" % sys.argv[0]
    sys.exit(0)

incontigfile = open(sys.argv[1])
distalPairs = open(sys.argv[2])
outpathfile = open(sys.argv[3],'w')
outcontigfile = open(sys.argv[4],'w')
pathPrefix= 'RNAPATH'
overlap = 30

if '-prefix' in sys.argv:
    pathPrefix = sys.argv[sys.argv.index('-prefix') + 1]

if '-overlap' in sys.argv:
    overlap = int(sys.argv[sys.argv.index('-overlap') + 1])
    
outheader = '#settings: %s' % ' '.join(sys.argv)
print outheader
outpathfile.write(outheader + '\n')
contigDict = {}
graphDict = {}
detailDict = {}
currentChrom = ''
seq = ''
origSize = []
origSize2 = []
nameList = []
contigNum = 0

#get the contigs
for line in incontigfile:
    if '>' in line:
        if currentChrom !='':
            nameList.append(currentChrom)
            contigDict[contigNum] = seq
            origSize.append(len(seq))
            origSize2.append(len(seq))
            contigNum += 1
        currentChrom = line.strip().split()[0][1:]
        seq = ''
    else:
        seq += line.strip()

# get the original N50 
totalSize = sum(origSize)
halfSize = totalSize / 2
origSize2.sort()
origSize2.reverse()
cumLength = 0
for asize in origSize2:
    if cumLength + asize > halfSize:
        print "#contigs", contigNum
        print "N50", asize
        break
    cumLength += asize
print origSize2[:50]

#detailDict = [] * contigNum
edgeArray = []
edgeSenseDict = {}

def visitLink(fromVertex, ignoreList = []):
    returnPath = [fromVertex]
    toVertex = []
    for toindex in xrange(contigNum):
        if edgeArray[fromVertex][toindex] > 1 and toindex not in ignoreList:
            toVertex.append(toindex)
    for vertex in toVertex:
        if sum(edgeArray[vertex]) == edgeArray[fromVertex][vertex]:
            edgeArray[fromVertex][vertex] = 0
            edgeArray[vertex][fromVertex] = 0
            return returnPath + [vertex]
        else:
            edgeArray[fromVertex][vertex] = 0
            try:
                return returnPath + visitLink(vertex, returnPath)
            except:
                return returnPath + [vertex]
    return ''
print "building the adjacency graph"
# build the edge graph
initList = [0] * contigNum

edgeArray = zeros((contigNum, contigNum), int16)
del initList

print len(edgeArray)
print len(edgeArray[50])

contigToRowLookup = {}
verticesWithEdges = {}
vertexEdges = {}
notSoloDict = {}
print "processing distal pairs"
for line in distalPairs:
    if line[0] == '#':
        continue
    fields = line.strip().split()
    contA = 'chr' + fields[1]
    try:
        contig1 = contigToRowLookup[contA]
    except:
        try:
            contig1 = nameList.index(contA)
            contigToRowLookup[contA] = contig1
        except:
            print "problem with end1: ", line
            continue
    start1 = int(fields[2])
    sense1 = fields[3]

    contB = 'chr' + fields[4]
    try:
        contig2 = contigToRowLookup[contB]
    except:
        try:
            contig2 = nameList.index(contB)
            contigToRowLookup[contB] = contig2
        except:
            print "problem with end2: ", line
            continue
    start2 = int(fields[5])
    sense2 = fields[6]
    
    edgeArray[contig1][contig2] += 1
    edgeArray[contig2][contig1] += 1
    verticesWithEdges[contig1] = ''
    verticesWithEdges[contig2] = ''
    if (contig1, contig2) in edgeSenseDict:
        edgeSenseDict[contig1, contig2].append((sense1, sense2))
    elif (contig2, contig1) in edgeSenseDict:
        edgeSenseDict[contig2, contig1].append((sense2, sense1))
    else:
        edgeSenseDict[contig1, contig2] = [(sense1, sense2)]
    if contig1 in vertexEdges:
        if contig2 not in vertexEdges[contig1]:
            vertexEdges[contig1].append(contig2)
    else:
        vertexEdges[contig1] = [contig2]
    if contig2 in vertexEdges:
        if contig1 not in vertexEdges[contig2]:
            vertexEdges[contig2].append(contig1)
    else:
        vertexEdges[contig2] = [contig1]
        
    if edgeArray[contig1][contig2] > 1:
        notSoloDict[contig1] = ''
        notSoloDict[contig2] = ''
	
    #detailDict[contig1].append((contig2, start1, start2))
    #detailDict[contig2].append((contig1, start2, start1))

zeroedEdge = 0
willVisitList = verticesWithEdges.keys()
willVisitList.sort()
print "visiting %d vertices" % len(willVisitList)

print "cleaning up graph of edges with weight 1"
for rindex in willVisitList:
    if rindex not in notSoloDict:
        cindex = vertexEdges[rindex][0]
        edgeArray[rindex][cindex] = 0
        edgeArray[cindex][rindex] = 0
        zeroedEdge += 1
        del verticesWithEdges[rindex]
print "%d 1-edges zeroed out" % zeroedEdge

zeroedEdge = 0
willVisitList = verticesWithEdges.keys()
willVisitList.sort()
print "visiting %d vertices" % len(willVisitList)

leafList = []
print "picking top 2 edges per vertex - zero out others"
for rindex in willVisitList:
    vertices = vertexEdges[rindex]
    rEdges = []
    for avertex in vertices:
        if avertex in willVisitList:
            rEdges.append((edgeArray[rindex][avertex], avertex))
    if len(rEdges) > 2:
        rEdges.sort()
        rEdges.reverse()
        #print rEdges
        zeroedEdge += len(rEdges[2:])
        for (weight, cindex) in rEdges[2:]:
            edgeArray[rindex][cindex] = 0
            edgeArray[cindex][rindex] = 0
    elif len(rEdges) == 1:
        if edgeArray[rindex][rEdges[0][1]] > 1:
            leafList.append(rindex)
print "zeroed out %d lower-weight edges at vertices with degree > 2" % zeroedEdge

pathList = []
visitedDict = {}
leafList.sort()
print "traveling through the graph"
# only start at leaves
for rindex in leafList:
    path = ''
    try:
        visitedDict[rindex]
        pass
    except:
        path = visitLink(rindex)
        try:
            if len(path) > 1:
                for vertex in path:
                    visitedDict[vertex] = ''
                print path
                pathList.append(path)
        except:
            print "problem with path from %d: %s" % (rindex, str(path))
            pass
print "found %d paths" % len(pathList)			

newSizeList = []
pathID = 0
for path in pathList:
    pathID += 1
    outpathfile.write('chr%s%d: %s\n' % (pathPrefix, pathID, str(path))) 
    pathDescription = ''
    for vertex in path:
        pathDescription += nameList[vertex] + ','
    pathDescription = pathDescription[:-1]
    outpathfile.write(pathDescription + '\n' )
    currentVertex = path[0]
    currentSense = '+'
    assemblyList = currentVertex
    sequence = contigDict[currentVertex]
    for nextVertex in path[1:]:
        if (currentVertex, nextVertex) in edgeSenseDict:
            senseList = edgeSenseDict[currentVertex, nextVertex]
            FR = senseList.count(('+','-'))
            RF = senseList.count(('-','+'))
        else:
            senseList = edgeSenseDict[nextVertex, currentVertex]
            # flip
            FR = senseList.count(('-','+'))
            RF = senseList.count(('+','-'))
        FF = senseList.count(('+','+'))
        RR = senseList.count(('-','-'))
        if currentSense == '-':
            # we had flipped the upstream piece! Must flip again
            temp1 = FR
            temp2 = FF
            FR = RR
            FF = RF
            RR = temp1
            RF = temp2
        if FR >= FF and FR >= RR and FR >= RF:
            # we have FR - leave alone
            sense1 = '+'
            sense2 = '-'
            assemblyList = ((assemblyList, '+'), (nextVertex, '+'))
            seqleft = sequence[-20:]
            seqright = contigDict[nextVertex][:overlap]
            if seqleft in seqright:
                pos = seqright.index(seqleft)
                offset = pos + 20
                outstring = 'stitching %d and %d using %d overlap' % (currentVertex, nextVertex, offset)
                print outstring
                outpathfile.write(outstring + '\n')
                sequence += contigDict[nextVertex][offset:]
            else:
                sequence += 'NN' + contigDict[nextVertex]
            currentSense = '+'
        elif FF >= RR and FF >= RF:
            # we have FF - flip seqright
            sense1 = '+'
            sense2 = '+'
            assemblyList = ((assemblyList, '+'), (nextVertex, '-'))
            seqleft = sequence[-20:]
            seqright = complement(contigDict[nextVertex])[:overlap]
            if seqleft in seqright:
                pos = seqright.index(seqleft)
                offset = pos + 20
                outstring = 'stitching %d and %d using %d overlap' % (nextVertex, currentVertex, offset)
                print outstring
                outpathfile.write(outstring + '\n')
                sequence += complement(contigDict[nextVertex])[offset:]
            else:
                sequence += 'NN' + complement(contigDict[nextVertex])
            currentSense = '-'
        elif RR >= RF:
            # we have RR - flip seqleft
            sense1 = '-'
            sense2 = '-'
            assemblyList = ((assemblyList, '-'), (nextVertex, '+'))
            seqleft = complement(sequence)[:20]
            seqright = contigDict[nextVertex][:overlap]
            if seqleft in seqright:
                pos = seqright.index(seqleft)
                offset = pos + 20
                outstring = 'stitching %d and %d using %d overlap' % (nextVertex, currentVertex, offset)
                print outstring
                outpathfile.write(outstring + '\n')
                sequence = complement(sequence) + contigDict[nextVertex][offset:]
            else:
                sequence = complement(sequence) + 'NN' + contigDict[nextVertex]
            currentSense = '+'
        else:
            # we have RF - flip both
            sense1 = '-'
            sense2 = '+'
            assemblyList = ((assemblyList, '-'), (nextVertex, '-'))
            seqleft = complement(sequence)[-20:]
            seqright = complement(contigDict[nextVertex])[:overlap]
            if seqleft in seqright:
                pos = seqright.index(seqleft)
                offset = pos + 20
                outstring = 'stitching %d and %d using %d overlap' % (nextVertex, currentVertex, offset)
                print outstring
                outpathfile.write(outstring + '\n')
                sequence = complement(sequence) + complement(contigDict[nextVertex])[offset:]
            else:
                sequence = complement(sequence) + 'NN' + complement(contigDict[nextVertex])
            currentSense = '-'
        outstring = '(%d, %d): FF %d RR %d RF %d FR %d : %s %s\t%s' % (currentVertex, nextVertex, FF, RR, RF, FR, sense1, sense2, str(assemblyList))
        print outstring
        outpathfile.write(outstring + '\n') 
        currentVertex = nextVertex
    outcontigfile.write('>chr%s%d %dbp %s | %s\n%s\n' % (pathPrefix, pathID, len(sequence), pathDescription, str(assemblyList), sequence))
    newSizeList.append(len(sequence))
            
for vertex in contigDict:
    if vertex in visitedDict:
        continue
    newSizeList.append(len(contigDict[vertex]))
    outcontigfile.write('>%s\n%s\n' % (nameList[vertex], contigDict[vertex]))

newSizeList.sort()
newSizeList.reverse()
cumLength = 0
for asize in newSizeList:
    if cumLength + asize > halfSize:
        print "#contigs", len(newSizeList)
        print "N50", asize
        break
    cumLength += asize
print newSizeList[:50]
