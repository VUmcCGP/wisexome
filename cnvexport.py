import numpy
import math

def writeToPickle(fileDrop, regional, tChrom, byRelative, byRegion, exonInfo, filteredPost, means, nonOccRegion):
	'Most important feature'
	import pickle
	exportList=[]
	regional.sort()
	for i,region in enumerate(regional):
		exportDict=dict()
		exportDict['chromosome']=tChrom
		exportDict['probes']=region[1]-region[0]+1
		exportDict['median']=numpy.median([x for x in byRelative[region[0]:region[-1]+1] if not math.isnan(x)])

		exportDict['minStart']=byRegion[region[0]][0]
		exportDict['minEnd']=byRegion[region[1]][1]

		exportDict['maxStart']=byRegion[region[0]-1][1]
		#exportDict['maxEnd']=byRegion[region[1]+1][0]
		exportDict['maxEnd']=byRegion[min(len(byRegion)-1,region[1]+1)][0]
		exportDict['filteredPost']=filteredPost[i]
		exportDict['mean']=means[i]
		exportDict['nonOcc']=nonOccRegion[i]

		hitGenes=[]
		for i,val in enumerate(exonInfo):
			if val[1] > exportDict['maxStart'] and val[0] < exportDict['maxEnd']:
				hitGenes.append(val[-1])
		exportDict['genes']=set(hitGenes)
		exportList.append(exportDict)

	pickle.dump(exportList,open(fileDrop+'.pickle','wb'))



def writeToBed(fileDrop):
	'Not finished'
	with open(fileDrop+'.bed', 'w') as toBedFile:
		lastProbe=0
		for i,val in enumerate(byRegion):
			bedScore=byRelative[i]#999
			#if byRelative[i] <= 999:
			#	bedScore=int((byRelative[i]+1)*2+0.5)
			bedScore = min(bedScore,999)

			if math.isnan(bedScore):
				bedScore=0

			colorBED=colorPalette[0]
			if i in extCall:
				#print bedScore

				if byRelative[i] > 0:
					colorBED=colorPalette[5]
				else:
					colorBED=colorPalette[6]

				if i in ignoreBins:
					if byRelative[i] > 0:
						colorBED=colorPalette[2]
					else:
						colorBED=colorPalette[1]

			else:
				colorBED=(0.5,0.5,0.5)#colorPalette[6]
				if i in ignoreBins:
					colorBED=(0.7,0.7,0.7)#colorPalette[1]


			bedProbePos=[]
			bedProbeLen=[]
			for probeIndex in range(lastProbe,len(probeInfo)):
				probeVal=probeInfo[probeIndex]
				if probeVal[0] > val[1]:
					break
				if probeVal[0] < val[1] and probeVal[1] > val[0]:
					bedProbePos.append(str(probeVal[0]-val[0]))
					bedProbeLen.append(str(probeVal[1]-probeVal[0]))
				lastProbe+=1

			colorRGB=[str(int(x*255.)) for x in colorBED]
			#print bedScore
			outLine=[tChrom,val[0],val[1],val[-1],bedScore,'+',val[0],val[0],','.join(colorRGB)]#,len(bedProbePos),','.join(bedProbeLen),','.join(bedProbePos)]
			toBedFile.write(' '.join([str(x) for x in outLine])+'\n')

def writeToBedShort(fileDrop):
	'Not finished'
	with open(fileDrop+'_short.bed', 'w') as toBedFile:
		for region in regional:
			posStart=byRegion[region[0]-1][1]
			print len(byRegion),region[1]+2,min(len(byRegion)-1,region[1]+2),byRegion[min(len(byRegion)-1,region[1]+2)][0]
			#posEnd=byRegion[region[1]+2][0]
			posEnd=byRegion[min(len(byRegion)-1,region[1]+2)][0]
			#print region,region[1]-region[0]+1,posStart,posEnd
				#print exon,byRegion[exon][0],byRegion[exon][1]
			exonLengths=['0']
			exonStarts=['0']

			for exon in range(region[0],region[1]+1):
				exonLengths.append(str(byRegion[exon][1]-byRegion[exon][0]))
				exonStarts.append(str(byRegion[exon][0]-posStart))


			exonLengths.append('0')
			exonStarts.append(str(posEnd-posStart))

			outLine=tChrom,posStart,posEnd,'NameTag',0,'+',byRegion[region[0]][0],byRegion[region[1]][1],'255,0,0',region[1]-region[0]+1,','.join(exonLengths),','.join(exonStarts)
				#byRegion[i][0], byRegion[i][1]]
			toBedFile.write(' '.join([str(x) for x in outLine])+'\n')


	#SampleID, Chr, Start, End, Value
def writeToSeg(fileDrop):
	'Not finished'
	SampleID=testData.split('/')[-1].split('.')[0]
	#tChrom

	with open(fileDrop+'.seg', 'w') as toSegFile:
		#toSegFile.write('SampleID\tChr\tStart\tEnd\tValue\n')
		for region in regional:
			noNan = [x for x in byRelative[region[0]:region[-1]+1] if not math.isnan(x)]
			if len(noNan) <= 2:
				continue
			segScore=numpy.median(noNan)
			posStart=byRegion[region[0]-1][1]
			#posEnd=byRegion[region[1]+2][0]
			posEnd=byRegion[min(len(byRegion)-1,region[1]+2)][0]

			segScore=int((segScore+1.25)*2)

			outLine=[SampleID,tChrom,posStart,posEnd,segScore]
			print '\t'.join([str(x) for x in outLine])
			toSegFile.write('\t'.join([str(x) for x in outLine])+'\n')
