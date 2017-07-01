from pylab import get_cmap
from matplotlib import gridspec
import matplotlib.pyplot as plt
import numpy
import math

colorPalette=[
	(0,0,0), # 0 black
	(0.90,0.60,0), # 1 Orange
	(0.35,0.70,0.90), # 2 Sky blue
	(0,0.60,0.50), # 3 Bluish green
	(0.95,0.90,0.25), # 4 Yellow
	(0,0.45,0.70), # 5 Blue
	(0.80,0.40,0), # 6 Vermillion
	(0.80,0.60,0.70), # 7 Reddish purple
	]

colorCal=colorPalette[5]#5
colorCal2=colorPalette[2]#2
colorFor=colorPalette[5]#6
colorFor2=colorPalette[2]#1
colorFor3=colorPalette[7]#1
colorNorm='gray'
colorNorm2='lightgray'
colorNaN=colorPalette[6]

colorHelp=(0.7,0.7,0.7)
colorThresh=colorPalette[3]#'green'

tendencyWindowSize=25

fColorMap='RdBu'#'seismic_r'


def drawZHelp(thresHold):
	plt.axhline(y=0, linewidth=0.75, color=colorHelp)
	plt.axhline(y=thresHold, linewidth=0.50, color=colorThresh)
	plt.axhline(y=-thresHold, linewidth=0.50, color=colorThresh)


def drawRHelp(relThresh):
	plt.axhline(y=0, linewidth=1, color=colorHelp)
	plt.axhline(y=relThresh, linewidth=0.750, color=colorThresh)
	plt.axhline(y=-relThresh, linewidth=0.750, color=colorThresh)


	plt.axhline(y=1, linewidth=1, color=colorHelp)
	plt.axhline(y=0.5, linewidth=0.5, color=colorHelp)
	plt.axhline(y=-0.5, linewidth=0.5, color=colorHelp)
	plt.axhline(y=-1, linewidth=1, color=colorHelp)

def flamePlot(mapZ, noNanZ, markRegions, flatten, colorCeil):
	print '\nPlotting flameview'

	ratios=[4,1,4]
	gs = gridspec.GridSpec(len(ratios), 1, height_ratios=ratios)
	plt.figure(figsize=(24,9))


	# Direct
	plt.subplot(gs[0])

	plt.imshow(mapZ, cmap=get_cmap(fColorMap), interpolation='nearest', aspect='auto')
	plt.clim(-colorCeil,colorCeil)

	plt.title('Overview')
	plt.xlabel('Probe number')
	plt.ylabel('Window Step')


	# Markers
	plt.subplot(gs[1])

	flatMarkers = [0]*len(markRegions[0])
	for i, val in enumerate(markRegions[0]):
		if abs(val) > flatten:
			flatMarkers[i] = val
	markRegions.append(flatMarkers)

	plt.imshow(markRegions, cmap=get_cmap(fColorMap), interpolation='nearest', aspect='auto')
	plt.clim(-colorCeil,colorCeil)

	plt.title('Maximum values')
	plt.xlabel('Probe number')
	plt.ylabel('Post|Pre')

	# Flattened
	for i in range(len(mapZ)):
		for j in range(len(noNanZ)):
			stouff=mapZ[i][j]
			if abs(stouff) < flatten:
				stouff = 0
			mapZ[i][j]=stouff

	plt.subplot(gs[2])
	plt.imshow(mapZ, cmap=get_cmap(fColorMap), interpolation='nearest', aspect='auto')
	plt.clim(-colorCeil,colorCeil)

	plt.title('Filtered at abs(z) > ' + str(flatten))
	plt.xlabel('Probe number')
	plt.ylabel('Window Step')

	plt.tight_layout()
	return plt



def overviewPlot(byTarget, refStdDev, byRelative, thresHold, relThresh):
	print '\nPlotting overview'

	def stoufferLine(data,window,colorChoice):
		tendency=[]
		for i in range(len(data)):
			theSet=data[max(0,i-window/2):min(len(data),i+window/2+1)]
			theSet.sort()
			theSet=theSet[window/10+1:-(window/10+1)]
			#theSet=[x for x in theSet if not math.isnan(x)]
			theValue=sum(theSet)/numpy.sqrt(len(theSet))
			tendency.append(theValue)
		plt.plot(tendency,'-',color=colorChoice,zorder=7)

	def makeTendencyLine(data,window,colorChoice):
		tendency=[]
		for i in range(len(data)):
			theSet=data[max(0,i-window/2):min(len(data),i+window/2+1)]
			theValue=numpy.median(theSet)
			tendency.append(theValue)
		plt.plot(tendency,'-',color=colorChoice,zorder=7)

	plt.figure(figsize=(24,8))

	plt.subplot(211)
	drawZHelp(thresHold)
	scale=30


	colorMap=[min(0.5,max(x,0)) for x in refStdDev]
	plt.scatter(range(len(byTarget)),byTarget, marker='x', c=colorMap, cmap=plt.cm.coolwarm)

	stoufferLine(byTarget,tendencyWindowSize,colorPalette[6])
	plt.ylim([-scale,scale])
	plt.xlim([-len(byTarget)*0.01,len(byTarget)+len(byTarget)*0.01])

	plt.title('Overview')
	plt.xlabel('Probe number')
	plt.ylabel('z-Score')
	plt.tight_layout()

	plt.subplot(212)
	drawRHelp(relThresh)
	scale=2

	for i,val in enumerate(byRelative):
		if not numpy.isfinite(val):
			#print i,val
			byRelative[i]=0

	plt.scatter(range(len(byRelative)),byRelative, marker='x', c=colorMap, cmap=plt.cm.coolwarm)

	makeTendencyLine(byRelative,tendencyWindowSize,colorPalette[1])
	plt.ylim([-scale,scale])
	plt.xlim([-len(byRelative)*0.01,len(byRelative)+len(byRelative)*0.01])

	plt.xlabel('Probe number')
	plt.ylabel('Copy number change')

	return plt



def regionPlot(tChrom,region,miniCall,calledIndex,mapZ,mapZMaxes,relative,byRelative,byRegion,ignoreBins,noNanC,noNanE,exonInfo,probeInfo,flatten,colorCeil,relThresh):
	ratios=[5,1,1,6,16]
	gs = gridspec.GridSpec(len(ratios), 1, height_ratios=ratios)
	plt.figure(figsize=(24,10))


	plt.subplot(gs[0])
	#miniCall=miniCalls[calledIndex]
	miniVals=[]
	length=miniCall[1]-miniCall[0]+1
	leftHang = min(miniCall[0],length)
	rightHang = min(len(mapZ[0])-1-miniCall[1],length+1)

	print leftHang,miniCall[0],miniCall[1]-miniCall[0]+1,miniCall[1],rightHang

	print miniCall[0]-length,miniCall[1]+length
	for i in range(len(mapZ)):
		miniVals.append(mapZ[i] [miniCall[0]-leftHang:miniCall[1]+rightHang])

	plt.imshow(miniVals,cmap=get_cmap(fColorMap), interpolation='nearest',aspect='auto')#interpolation='nearest',
	plt.clim(-colorCeil+flatten,colorCeil-flatten)
	#plt.colorbar(orientation='vertical')
	plt.tight_layout()

	plt.subplot(gs[1])
	plt.imshow([[x[1] for x in mapZMaxes[miniCall[0]-leftHang:miniCall[1]+rightHang]]],cmap=get_cmap(fColorMap), interpolation='nearest',aspect='auto')#interpolation='nearest',
	plt.clim(-colorCeil+flatten,colorCeil-flatten)
	#plt.colorbar(orientation='vertical')
	plt.tight_layout()

	plt.subplot(gs[2])

	#marker=[0]*min(miniCall[0],length) + [1]*(miniCall[1]-miniCall[0]) + [0]*min(len(mapZ[0])-1-miniCall[1],length)
	marker=[0]*leftHang + [1]*(miniCall[1]-miniCall[0]+1) + [0]*(rightHang-1)


	plt.imshow([marker],cmap=get_cmap(fColorMap), interpolation='nearest',aspect='auto')#interpolation='nearest',
	plt.clim(-1,1)
	#plt.colorbar(orientation='vertical')
	plt.tight_layout()

	plt.subplot(gs[3])
	plt.plot(noNanC[miniCall[0]-leftHang:miniCall[1]+rightHang])
	plt.plot(noNanE[miniCall[0]-leftHang:miniCall[1]+rightHang])
	#plt.xlim([miniCall[0]-leftHang,miniCall[1]+rightHang-1])
	plt.xlim([0,len(noNanC[miniCall[0]-leftHang:miniCall[1]+rightHang])-1])

	plt.subplot(gs[-1])
	noNan = [x for x in byRelative[region[0]:region[-1]+1] if not math.isnan(x)]
	
	'''
	roundedEffect=numpy.median(byRelative[region[0]:region[-1]+1])
	if roundedEffect>0:
		roundedEffect+=0.25
	else:
		roundedEffect-=0.25
	roundedEffect=int(roundedEffect*2)/2.0
	print roundedEffect,numpy.median(noNan)
	'''
	scale=1.25

	#drawRHelp()
	drawRHelp(relThresh)

	regStart=byRegion[region[0]][0]
	regEnd=byRegion[region[-1]][1]
	regWidth=regEnd-regStart
	regHang=max(0.3*regWidth,10000)
	plt.xlim([regStart-regHang, regEnd+regHang])
	print 'Plotting region:',regStart-regHang, regEnd+regHang
	# Surrounding non called probes
	# Connecting lines
	surroundColor=colorNorm
	for i in range(1,len(byRegion)):
		if byRegion[i][0] > regStart-regHang and byRegion[i-1][1] < regEnd+regHang:
			plt.plot([byRegion[i-1][1], byRegion[i][0]], [byRelative[i-1],byRelative[i]], color=surroundColor, linestyle='-', linewidth=1)

	# Add missing lines due to NaNs
	lastVal=byRelative[0]
	lastPos=byRegion[0][1]
	for i in range(1,len(byRegion)):
		if byRegion[i][0] > regStart-regHang and byRegion[i-1][1] < regEnd+regHang:
			if math.isnan(byRelative[i-1]) and not math.isnan(byRelative[i]):
				plt.plot([lastPos, byRegion[i][0]], [lastVal,byRelative[i]], color=surroundColor, linestyle='--', linewidth=1)
		if not math.isnan(byRelative[i]):
			lastVal=byRelative[i]
			lastPos=byRegion[i][1]

	# The probes
	for i in range(0,len(byRegion)):
		if byRegion[i][1] > regStart-regHang and byRegion[i][0] < regEnd+regHang:
			if i not in ignoreBins:
				plt.plot([byRegion[i][0], byRegion[i][1]], [byRelative[i],byRelative[i]], color=surroundColor, linestyle='-', linewidth=4)
			else:
				plt.plot([byRegion[i][0], byRegion[i][1]], [byRelative[i],byRelative[i]], color=colorNorm2, linestyle='-', linewidth=4)

	# Probes that are NaNs
	for i in range(0,len(byRegion)):
		if byRegion[i][1] > regStart-regHang and byRegion[i][0] < regEnd+regHang:
			if math.isnan(byRelative[i]):
				plt.plot([byRegion[i][0], byRegion[i][1]], [0,0], color=colorNaN, linestyle='-', linewidth=4)

	# Called probes
	calledColor=colorCal
	calledOccColor=colorCal2

	regionMedian = numpy.median(byRelative[region[0]:region[1]+1])
	surroundWidth=6

	# Actual exon position marking and annotation
	prevAnn=''
	annCount=0
	probeFinger=0
	annStart=0
	annEnd=0

	for i,val in enumerate(exonInfo):
		newAnn = val[-1]
		if val[1] > regStart-regHang and val[0] < regEnd+regHang:

			shape=[]

			while probeInfo[probeFinger][1] > val[0] and probeFinger > 0:
				probeFinger -= 1

			while probeInfo[probeFinger][1] < val[0]  and probeFinger < len(probeInfo)-1:
				probeFinger += 1

			while probeInfo[probeFinger][0] < val[1]  and probeFinger < len(probeInfo)-1:
				thisProbe=[probeInfo[probeFinger][0],probeInfo[probeFinger][1],relative[probeFinger]]

				shape.append(thisProbe)
				probeFinger += 1

			if val[1] > regStart and val[0] < regEnd:
				plt.plot(val[:2], [0,0], color=colorPalette[1], linestyle='-', linewidth=4,zorder=4)#[roundedEffect,roundedEffect]
			else:
				plt.plot(val[:2], [0,0], color=colorPalette[6], linestyle='-', linewidth=4,zorder=3)

			if newAnn != prevAnn:
				plt.plot([annStart,annEnd], [0,0], color=colorPalette[0], linestyle='-', linewidth=3,zorder=2)

				annStart = val[0]
				upside=1
				offSet=0
				if sum(byRelative[region[0]:region[1]+1])>0:
					upside=-1
				if annCount%2==1:
					offSet=0.1
				plt.annotate(newAnn, xy=(val[0], 0), xytext=(val[0], +(0.75+offSet)*upside-0.05),#byRelative[i]+0.05))
					arrowprops=dict(facecolor='lightgray', shrink=0.1))
				annCount+=1
			annEnd = val[1]
			prevAnn = newAnn

	plt.plot([annStart,annEnd], [0,0], color=colorPalette[0], linestyle='-', linewidth=3,zorder=2)

	# Connecting lines
	for i in range(region[0]+1,region[1]+1):
		if math.isnan(byRelative[i-1]) \
		or byRelative[i-1] == 0.0 \
		or math.isnan(byRelative[i]) \
		or byRelative[i] == 0.0 \
		or i in ignoreBins \
		or i-1 in ignoreBins:
			plt.plot([byRegion[i-1][1], byRegion[i][0]], [byRelative[i-1],byRelative[i]], color=calledOccColor, linestyle='--', linewidth=2)
		else:
			plt.plot([byRegion[i-1][1], byRegion[i][0]], [byRelative[i-1],byRelative[i]], color=calledColor, linestyle='-', linewidth=2)

	# Add missing lines due to NaNs
	lastVal=byRelative[region[0]]
	lastPos=byRegion[region[0]][1]
	for i in range(region[0]+1,region[1]+1):
		if (math.isnan(byRelative[i-1]) or byRelative[i-1] == 0.0 or i-1 in ignoreBins) \
		and not (math.isnan(byRelative[i]) or byRelative[i] == 0.0 or i in ignoreBins):
			plt.plot([lastPos, byRegion[i][0]], [lastVal,byRelative[i]], color=calledColor, linestyle='-', linewidth=2)
		if not (math.isnan(byRelative[i]) or byRelative[i] == 0.0 or i in ignoreBins):
			lastVal=byRelative[i]
			lastPos=byRegion[i][1]

	# The probes
	for i in range(region[0],region[1]+1):
		if i not in ignoreBins:
			if math.isnan(byRelative[i]) or byRelative[i] == 0.0:
				plt.plot([byRegion[i][0], byRegion[i][1]], [byRelative[i],byRelative[i]], color=calledOccColor, linestyle='-', linewidth=6)
			else:
				plt.plot([byRegion[i][0], byRegion[i][1]], [byRelative[i],byRelative[i]], color=calledColor, linestyle='-', linewidth=6)
		else:
			plt.plot([byRegion[i][0], byRegion[i][1]], [byRelative[i],byRelative[i]], color=calledOccColor, linestyle='-', linewidth=6)


	plt.title('CNV: chr'+str(tChrom)+':'+str(regStart)+'-'+str(regEnd)+' ('+str((regEnd-regStart)/1000)+' Kbp)')
	#if markRegion != '' and region == regional[-1]:
	#	plt.title('Marked: '+str(tChrom)+':'+str(regStart)+'-'+str(regEnd)+' ('+str((regEnd-regStart)/1000)+' Kbp)')

	plt.xlabel('bp Position on chromosome')
	plt.ylabel('Copy number change')

	maxCNV=max([byRelative[x] for x in range(region[0],region[-1]+1) if x not in ignoreBins])

	plt.ylim([-1.25,max(1.25,maxCNV*1.25)])

	return plt

	'''
		upDown='dup'
		if regionMedian < 0:
			upDown = 'del'

		if makePlotFile:
			print 'Writing plots to file'
			plt.savefig(fileDrop+'_'+str(regStart)+'-'+str(regEnd)+'.'+upDown+'.pdf')

	if makePlot:
		print 'Showing plots'
		plt.show()
	'''
