import cPickle
import argparse
import numpy
import warnings
import math
import matplotlib

#import scipy
#import scipy.stats


# Some default variable values
thresHold = 5.64
relThresh = 0.35
connectZMax = 1
maxRange = 8#16#8


def getZScore(value, reference):
	average = numpy.average(reference)
	stddev = numpy.std(reference)

	if stddev == 0:
		return 0

	Z = (value - average) / stddev

	return Z


def getReference(lookUp, cutOff):
	reference = []
	removed = 0
	for exon in lookUp:
		# print exon
		if len(exon) > 0:
			if float(exon[0][0]) < cutOff:
				reference.append(float(exon[0][0]))
			else:
				removed += 1
	return reference, removed


def getOptimalCutoff(lookUp, repeats, optimalCutoff):
	for i in range(0, repeats):
		reference, removed = getReference(lookUp, optimalCutoff)
		average = numpy.average(reference)
		stddev = numpy.std(reference)
		optimalCutoff = average + 3 * stddev
	return optimalCutoff


# Data loaders

def loadOccurrences(fileOcc):
	print '\tLoading:\t' + fileOcc
	ignoreBins = []
	with open(fileOcc, 'r') as probeData:
		for line in probeData:
			splitLine = line.split()
			if splitLine[0] != 'Loading:' and splitLine[0] != '[]':
				if int(splitLine[1]) > 4:
					ignoreBins.append(int(splitLine[0]))
	return ignoreBins


def loadProbes(probeBed,tChrom):
	print '\tLoading:\t' + probeBed
	probeInfo = []
	with open(probeBed, 'r') as probeData:
		for line in probeData:
			splitLine = line.split()
			if splitLine[0][3:] == tChrom:
				start = int(splitLine[1])
				end = int(splitLine[2])
				probeName = splitLine[-1]
				probeInfo.append([start, end, probeName])
	probeInfo.sort()
	return probeInfo


def loadExons(exonBed,tChrom):
	print '\tLoading:\t' + exonBed
	exonInfo = []
	with open(exonBed, 'r') as exonData:
		for line in exonData:
			splitLine = line.split()
			if splitLine[2][3:] == tChrom and len(splitLine) > 10:
				exonCount = int(splitLine[8])
				if exonCount > 0:
					exonStarts = [int(x) for x in splitLine[9].split(',')[:-1]]
					exonEnds = [int(x) for x in splitLine[10].split(',')[:-1]]
					geneName = splitLine[12]
					for exonIndex in range(int(exonCount)):
						start = exonStarts[exonIndex]
						end = exonEnds[exonIndex]
						exonInfo.append([start, end, exonIndex, geneName])
	exonInfo.sort()
	return exonInfo


def loadSample(testData, tChrom):
	print '\tLoading:\t' + testData
	testSample = cPickle.load(open(testData, 'r'))

	if tChrom == 'X':
		sumCount = float(sum([sum(testSample[x]) for x in
							  testSample.keys()]))  # -sum(hitFile['chrX']) # - X is newly added after ex12
		targets = [x / sumCount for x in testSample['chr' + tChrom]]

		print [x for x in testSample['chrX'] if x > 0][:10]
		print numpy.median(targets)

		# Par regions: 60001-2699520; 154931044-155270560

		if numpy.median(targets) < 0.000002:
			print 'Patient seems Male'
			testSample['chrX'] = [x * 2 for x in testSample['chrX']]
		else:
			print 'Patient seems Female'

		print [x for x in testSample['chrX'] if x > 0][:10]
	return testSample


def loadReference(refData):
	print '\tLoading:\t' + refData
	reference = cPickle.load(open(refData, 'r'))
	return reference

def loadFilterBed(filtPostSoft,tChrom):
	print '\tLoading:\t' + filtPostSoft
	filtData = []
	with open(filtPostSoft, 'r') as filtFile:
		for line in filtFile:
			splitLine = line.split()
			if splitLine[0][3:] == tChrom:
				filtData.append([int(x) for x in splitLine[1:]])
	return filtData

# Data testing

def cnvTest(ignoreBins, probeInfo, testSample, reference, tChrom, refCutOff, directCalls):
	print 'Testing for CNVs'
	curChrom = testSample['chr' + str(tChrom)]
	results = []
	relative = []
	refSize = []
	refStdDev = []
	refMean = []

	print refCutOff
	for tarI, tarVal in enumerate(curChrom):
		refPlaces = reference[tarI]
		refSet = []
		for refI, refVal in enumerate(refPlaces):
			if refVal[0] < refCutOff:
				refChr = refVal[2]
				refPos = refVal[1]
				refSet.append(testSample['chr' + str(refChr)][refPos])
		if len(refSet) < 10:
			refSet = []
		relative.append(tarVal / numpy.mean(refSet))
		results.append(getZScore(tarVal, refSet))
		refSize.append(len(refSet))
		refStdDev.append(numpy.std(refSet) / (numpy.mean(refSet)))
		refMean.append(numpy.mean(refSet))

	print 'Refsizes min/max', min(refSize), max(refSize)

	byTarget = []
	byRelative = []
	curTarget = []
	curRelative = []
	lastProbeName = 'derp;derp'
	start = 0
	end = 0
	called = []
	byRegion = []

	# Probe direct based
	byTarget = results
	byRelative = [x - 1 for x in relative]
	byRegion = probeInfo

	for i in range(len(relative)):
		if abs(byTarget[i]) > thresHold and abs(byRelative[i]) > relThresh and i not in ignoreBins > 0:
			called.append(i)

	if directCalls:
		for call in called:
			print call, 'directCallTag'
		exit()

	noNanZ = []
	noNanR = []
	noNanI = []
	noNanC = []
	noNanE = []
	for i in range(len(byRelative)):
		if refMean[i] <= 10 or math.isnan(byRelative[i]) or i in ignoreBins:
			continue
		else:
			noNanZ.append(byTarget[i])
			noNanR.append(byRelative[i])
			noNanI.append(i)
			noNanC.append(curChrom[i])
			noNanE.append(refMean[i])

	print len(byTarget), len(noNanZ)  # ,len(noNanR),len(noNanI)

	mapZ = []  # [noNanZ]
	mapR = []  # [noNanR]
	for i in range(maxRange):
		print i, 'Working on window size:', i * 2 + 1
		tmpZ = [0] * len(noNanZ)
		tmpR = [0] * len(noNanZ)
		for j in range(len(noNanZ) - 1):
			leftEnd = max(0, j - i)
			rightEnd = min(len(noNanZ) - 1, j + i + 1)
			localData = noNanZ[leftEnd:rightEnd]
			stouff = sum(localData) / math.sqrt(len(localData))

			tmpZ[j] = stouff
			med = numpy.median(noNanR[leftEnd:rightEnd])
			tmpR[j] = (med)

		mapZ.append(tmpZ)
		mapR.append(tmpR)

	import operator
	mapZMaxes = [[] for x in range(len(mapZ[0]))]
	mapZT = map(list, zip(*mapZ))
	for i, val in enumerate(mapZT):
		maxIndex, maxVal = max(enumerate(val), key=operator.itemgetter(1))
		minIndex, minVal = min(enumerate(val), key=operator.itemgetter(1))
		if abs(maxVal) > abs(minVal):
			mapZMaxes[i] = [maxIndex, maxVal]
		else:
			mapZMaxes[i] = [minIndex, minVal]

	print sum([abs(x[1]) > thresHold for x in mapZMaxes])

	mapZMaxCalls = []
	for i, val in enumerate(mapZMaxes):
		if abs(val[1]) > thresHold:
			if False:
				continue
			else:
				mapZMaxCalls.append(i)

	mapZRegions = []
	if len(mapZMaxCalls) > 0:
		curStart = mapZMaxCalls[0]
		for i in range(1, len(mapZMaxCalls)):
			if mapZMaxCalls[i] - mapZMaxCalls[i - 1] > connectZMax:
				mapZRegions.append([curStart, mapZMaxCalls[i - 1]])
				curStart = mapZMaxCalls[i]
		mapZRegions.append([curStart, mapZMaxCalls[-1]])

	miniCalls = []
	for i, region in enumerate(mapZRegions):
		zMean=numpy.mean([x[1] for x in mapZMaxes[region[0]:region[1]+1]])
		print 'Current Region:',region,zMean
		outSide=[]

		# First everything to the right so positive values match right end data
		for j,otherRegion in enumerate(mapZRegions[i+1:]):
		#	print 'derpRight:',mapZRegions[i+j][1]+1,otherRegion[0]
			outSide.extend(noNanR[mapZRegions[i+j][1]+1:otherRegion[0]])

		# And the very end
		#print 'derpRight:',mapZRegions[-1][-1]+1,':'
		outSide.extend(noNanR[mapZRegions[-1][-1]+1:])

		# Then add everything to the left so negative values match left end data
		# Starting at zero
		#print 'derpLeft:',':',mapZRegions[0][0]
		outSide.extend(noNanR[:mapZRegions[0][0]])

		for j,otherRegion in enumerate(mapZRegions[:i]):
		#	print 'derpLeft:',otherRegion[1]+1,mapZRegions[j+1][0]
			outSide.extend(noNanR[otherRegion[1]+1:mapZRegions[j+1][0]])

		#print len(outSide),outSide[-100]

		inSide=noNanR[region[0]:region[1]+1]
		
		outSide = [min(1.,x) for x in outSide]
		inSide  = [min(1.,x) for x in inSide]

		reducedOutside=outSide[maxRange:-maxRange]
		#backEnd=
		#backEnd.reverse()
		extendedRegion=outSide[-maxRange:]+inSide+outSide[:maxRange]
		#print len(inSide),reducedOutside[0],reducedOutside[-1]
		#print len(extendedRegion),extendedRegion[0],extendedRegion[-1]

		#print "Lengths:",len(reducedOutside),len(extendedRegion)

		maxSegmentation = [0, 0, 0, 0]
		print 'Regionsize:',len(extendedRegion[maxRange:-maxRange])
		if region[1]-region[0]==0 and abs(inSide[0])>relThresh:
			maxSegmentation = [0.0000001,numpy.mean(inSide),region[0],region[1]]
		#print maxSegmentation
		
		#startPoint = max(region[0] - maxRange, 0)
		#endPoint = min(region[1] + maxRange + 1, len(noNanR) - 1)
		#print extendedRegion
		for x in range(0, len(extendedRegion) - maxRange):
			for y in range(max(x+1,maxRange), len(extendedRegion)):
				#left = numpy.mean(noNanR[startLeft:x])
				out = numpy.mean(reducedOutside+extendedRegion[:x]+extendedRegion[y:])
				mid = numpy.mean(extendedRegion[x:y])
				#right = numpy.mean(noNanR[y:startRight + 1])
				#diff = abs(left - mid) + abs(mid - right)
				diff = abs(mid-out)
				diff *= numpy.sqrt(y - x + 1)
				#print x,y,mid,out,diff
				'''
				meanA=out
				stdA=scipy.stats.tstd(reducedOutside+extendedRegion[:x]+extendedRegion[y:])
				obsA=len(reducedOutside+extendedRegion[:x]+extendedRegion[y:])
				meanB=mid
				#stdB=scipy.stats.tstd(extendedRegion[x:y])
				obsB=len(extendedRegion[x:y])
				
				diff = 1-scipy.stats.ttest_ind_from_stats(meanA,stdA,obsA,meanB,stdA,obsB)[1]
				'''
				#diff = 1-scipy.stats.ttest_ind(reducedOutside+extendedRegion[:x]+extendedRegion[y:],extendedRegion[x:y])[1]
				#print len(noNanZ[x+region[0]-maxRange : y+region[0]-maxRange - 1])
				#diff=sum(noNanZ[x+region[0]-maxRange : y+region[0]-maxRange - 1])/math.sqrt(y-x+1)

				# Region is not the right 'direction'
				if (mid > 0) != (zMean > 0):
					#print 'This does not make much sense:',mid,zMean
					diff = 0
				# Region does not meet our thresholds
				if abs(mid) < relThresh or abs(numpy.median(extendedRegion[x:y])) < relThresh:
					diff = 0
				# Don't bother with some end points if they don't fulfill our threshold anyway
				if abs(extendedRegion[x]) < relThresh or abs(extendedRegion[y-1]) < relThresh or \
					(extendedRegion[x] > 0) != (mid > 0) or (extendedRegion[y-1] > 0) != (mid > 0) :
					diff = 0 
				# Update our champion
				if diff > maxSegmentation[0]:
					maxSegmentation = [diff, mid, x+region[0]-maxRange, y+region[0]-maxRange - 1]
					
				#if diff > 0 and len(extendedRegion[x:y]) >= 2:
				#	print i,x,y,diff,scipy.stats.ttest_ind(reducedOutside+extendedRegion[:x]+extendedRegion[y:],extendedRegion[x:y])
					#print "Lengthseg:",len(reducedOutside+extendedRegion[:x]+extendedRegion[y:]),len(extendedRegion[x:y])
					#print extendedRegion[x:y]
		#print 'Segmentation:', maxSegmentation[-1] - maxSegmentation[-2] + 1, maxSegmentation
		#print 'Put back:',noNanI[maxSegmentation[-2]], noNanI[maxSegmentation[-1]]
		#print relative[noNanI[maxSegmentation[-2]-1] : noNanI[maxSegmentation[-1]+1]]
		#print noNanR[maxSegmentation[-2]-1 : maxSegmentation[-1]+1]
		#print '\n\n'

		if maxSegmentation[-1] - maxSegmentation[-2] > -1 and maxSegmentation[-2] > 0:
			miniCalls.append([maxSegmentation[-2], maxSegmentation[-1]])
			print maxSegmentation
			
	combinedCalls=[]
	popped=[]
	print miniCalls
	for i,call in enumerate(miniCalls):
		if i in popped:
			continue
		curCall=call[:]
		#print i,call
		#print popped
		for j,other in enumerate(miniCalls[i+1:]):
			#print '',j,other
			if other[0] <= curCall[1]:
				curCall[1] = other[1]
				popped.append(i+j+1)
				#print 'pop!'
		combinedCalls.append(curCall)
	print combinedCalls
	
	miniCalls = combinedCalls
	if miniCalls != []  and miniCalls[-1][1] >= len(noNanI):
		print "This call is fishy:",miniCalls[-1],"It ends beyond the last testable probes?"
		miniCalls[-1][1] = len(noNanI)-1
		
	for call in miniCalls:
		data=noNanZ[call[0]:call[1]+1]
		stouff = sum(data) / math.sqrt(len(data))
		print call,stouff
	regional = [[noNanI[x[0]], noNanI[x[1]]] for x in miniCalls]  # []
	means = [numpy.mean(noNanR[x[0]:min(len(noNanI)-1,x[1]+1)]) for x in miniCalls]

	#extCall = called[:]
	maxDistances = [noNanI[min(len(noNanI)-1,x[1]+1)]-noNanI[x[0]-1] for x in miniCalls]
	distances = [x[1]-x[0] for x in miniCalls]
	#print miniCalls,regional
	#print maxDistances,distances
	nonOccRegion=[(maxDistances[x]-distances[x]-2) for x in range(len(distances))]
	#print nonOccRegion

	return called,regional,miniCalls,mapZ,mapR,noNanZ,noNanR,mapZMaxes,byTarget, refStdDev, byRelative, relative, byRegion,noNanC,noNanE,nonOccRegion,means

def filterPostSoft(probeInfo,regional,filtPSRegions):
	filtered=[1]*len(regional)
	#print 'derp',filtPSRegions
	for i,region in enumerate(regional):
		callLeft = probeInfo[region[0]][0]
		callRight = probeInfo[region[1]+1][1]
		#print callLeft,callRight

		for j,target in enumerate(filtPSRegions):
			filtLeft=target[0]
			filtRight=target[1]

			#print callLeft,callRight,target

			if filtLeft <= callRight and filtRight >= callLeft:
				filtered[i] = 0
	return filtered

def main():
	parser = argparse.ArgumentParser(description='TBA (or CBA perhaps)',
									 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('testfile', type=str,
						help='file containing amount of hits per probe')
	parser.add_argument('reffile', type=str,
						help='file containing reference probes per probe target')
	parser.add_argument('dropfile', type=str,
						help='file to save output')
	parser.add_argument('targetchr', type=str,
						help='chromosome to target bins on')
	parser.add_argument('probefile', type=str,
						help='bed file with probe locations')
	parser.add_argument('exonfile', type=str,
						help='bed file with exon information')
	parser.add_argument('occfile', type=str,
						help='empty file or file with amount of calls per exon over previous runs')
	parser.add_argument('-filtpostsoft', type=str, default='',
						help='ignore any call not touching a target region in this bed file')
	parser.add_argument('-plot', action='store_true',
						help='plot the data interactively')
	parser.add_argument('-plotfile', action='store_true',
						help='plot the data to a file')
	parser.add_argument('-mark', type=str, default='',
						help='use start-end, no commas, used to force plot an area (blue)')
	parser.add_argument('-mpluse', default='agg', type=str,
						help='make matplotlib use another backend for plotting')
	parser.add_argument('-direct', action='store_true',
						help='use for training step to create occurrence file')

	args = parser.parse_args()

	testData = args.testfile
	refData = args.reffile
	probeBed = args.probefile
	exonBed = args.exonfile
	tChrom = args.targetchr
	fileDrop = args.dropfile
	fileOcc = args.occfile
	makePlot = args.plot
	makePlotFile = args.plotfile
	markRegion = args.mark
	directCalls = args.direct
	filtPostSoft = args.filtpostsoft

	matplotlib.use(args.mpluse)

	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		ignoreBins = loadOccurrences(fileOcc)
		probeInfo = loadProbes(probeBed, tChrom)
		exonInfo = loadExons(exonBed, tChrom)
		testSample = loadSample(testData, tChrom)
		reference = loadReference(refData)

		print 'Testing for reference cutoff value'
		refCutOff = getOptimalCutoff(reference, 3, 1)

		called,regional,miniCalls,mapZ,mapR,noNanZ,noNanR,mapZMaxes,byTarget, refStdDev, byRelative, relative, byRegion,noNanC,noNanE,nonOccRegion,means = cnvTest(ignoreBins, probeInfo, testSample, reference, tChrom, refCutOff, directCalls)

	#print called,regional,miniCalls

	markRegions = [[x[1] for x in mapZMaxes]]

	filteredPost = [0]*len(regional)
	if filtPostSoft is not '':
		filtPSRegions = loadFilterBed(filtPostSoft,tChrom)
		filteredPost = filterPostSoft(probeInfo,regional,filtPSRegions)
	#print regional
	#print miniCalls

	if makePlot or makePlotFile:
		import matplotlib.pyplot as plt
		import cnvplot
		flamePlot = cnvplot.flamePlot(mapZ, noNanZ, markRegions, thresHold, 15)
		overviewPlot = cnvplot.overviewPlot(byTarget, refStdDev, byRelative, thresHold, relThresh)

		regionPlots=[]
		for calledIndex,region in enumerate(regional):
			if filteredPost[calledIndex] == 1:
				continue

			#print ''
			miniCall = miniCalls[calledIndex]
			print calledIndex,region,nonOccRegion[calledIndex],
			if region[1] - region[0] - nonOccRegion[calledIndex] < 3: # or nonOccRegion[calledIndex] < 0.75:
				print 'skipped'
				#continue
			print 'kept'
			regionPlots.append([region[0],region[1],cnvplot.regionPlot(tChrom,region,miniCall,calledIndex,mapZ,mapZMaxes,relative,byRelative,byRegion,ignoreBins,noNanC,noNanE,exonInfo,probeInfo,thresHold, 15,relThresh)])
			print ''

		if makePlot:
			plt.show()
		if makePlotFile:
			plt.figure(1)
			plt.savefig(fileDrop+'.flameview.pdf')

			plt.figure(2)
			plt.savefig(fileDrop+'.overview.pdf')

			#print 'location',str(probeInfo[region[0][0]]),str(probeInfo[region[1][1]])

			for i,reg in enumerate(regionPlots):
				plt.figure(i+3) # Due to some undefined logic apparently we can't use figs added to the list previously
				plt.savefig(fileDrop+'_'+str(probeInfo[reg[0]][0])+'-'+str(probeInfo[reg[1]][1])+'.pdf')

	import cnvexport as exp

	exp.writeToPickle(fileDrop, regional, tChrom, byRelative, byRegion, exonInfo, filteredPost, means, nonOccRegion)


if __name__ == "__main__":
	main()
