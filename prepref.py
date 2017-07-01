import cPickle
import argparse
import glob
import time
import numpy

def quickSelect(bigList,amount):
	bestList=[[1]]*amount
	tmpList=[[1]]*amount
	curV = 0
	worst = 1
	for iR,valR in enumerate(bigList):
		if curV == amount:
			curV=0
			bestList = sorted(bestList + tmpList)[:amount]
			worst=max([x[0] for x in bestList])
		if valR < worst:
			tmpList[curV]=[valR,iR]
			curV+=1
	return sorted(bestList + tmpList[:curV])[:amount]
	
def notQuickSelect(bigList,amount):
	tmpList = [[val,i] for i,val in enumerate(bigList)]
	return sorted(tmpList)[:amount]

parser = argparse.ArgumentParser(description='TBA (or CBA perhaps)',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('refdir', type=str,
					help='files containg read hits per probe')
parser.add_argument('dropfile', type=str,
					help='file containg not too different probes')
parser.add_argument('targetchr', type=str,
					help='chromosome to target bins on')
parser.add_argument('refchr', type=str,
					help='chromosome to find reference bins on')
args = parser.parse_args()

refData		= args.refdir
fileDrop	= args.dropfile
tChr		= args.targetchr
rChr		= args.refchr

#sampleHits = cPickle.load(open(fileHits,'r'))

distances	= []

targets		= []
references	= []

referenceFiles = glob.glob(refData + '/*.hits')
for refFile in referenceFiles:
	print '\tLoading:\t' + refFile
	hitFile =  cPickle.load(open(refFile,'r'))
	#print hitFile.keys()
	sumCount = float(sum([sum(hitFile[x]) for x in hitFile.keys()]))#-sum(hitFile['chrX']) # - X is newly added after ex12
	target=[x/sumCount for x in hitFile[tChr]]
	reference=[x/sumCount for x in hitFile[rChr]]
	
	#targets.append(target)
	#references.append(reference)
	
	if tChr == 'chrX' and numpy.median(target) < 0.00000175:#0.000002:
		target=[x/sumCount*2 for x in hitFile[tChr]]
		# Well nevermind but technically we could accept the male samples for training
	else:
		targets.append(target)
		references.append(reference)
		
		
	#print "tagthis:",refFile,numpy.median(targets[-1]),numpy.median(references[-1]),sumCount
#	if tChr=='chrX':
#		mover=numpy.median(references)/numpy.median(targets)
#		targets = [x*mover for x in targets]

	#print "tagthis:",refFile,numpy.median(targets),numpy.median(references),sumCount

#exit()
#for x in targets[0]:
#	distances.append([0] * len(references[0]))

#print len(distances),'vs',len(distances[0])

tarsT = map(list, zip(*targets))
refsT = map(list, zip(*references))

print len(tarsT),len(tarsT[0])
print len(refsT),len(refsT[0])

t = time.time()
toDoLen=len(tarsT)
distances = [[]] * len(tarsT)
#diffLists = [[]] * len(tarsT)
for tarTi,tarTv in enumerate(tarsT):
	curDists=[0]*len(refsT)
	for refTi,refTv in enumerate(refsT):
		curDist=0
		for i in range(len(tarTv)):
			diff = tarTv[i]-refTv[i]
			curDist+=diff*diff
		curDists[refTi]=curDist
	distances[tarTi]=quickSelect(curDists,100)

	#curDiffs=[0]*len(distances[tarTi])
	#for selPi,selPval in enumerate(distances[tarTi]):
	#	curRef = refsT[selPval[1]]
	#	curSubDiff = [[]] * len(tarTv)
	#	for i in range(len(tarTv)):
	#		diff=tarTv[i]-curRef[i]
	#		curSubDiff[i]=diff*diff
	#	curDiffs[selPi]=curSubDiff
	#diffLists[tarTi]=curDiffs

	if tarTi % 50 == 0:
		elapsed = time.time() - t
		totalTime = elapsed/(float(tarTi+1)/toDoLen)
		print 'TargetProbe:',tarTi,'\tProgress:',int(float(tarTi+1)/toDoLen*100),'%\tTimeElapsed:',int(elapsed)/60,'m',int(elapsed)%60,'s\tTimeLeft:',int((totalTime-elapsed)/60)/60,'h',int((totalTime-elapsed)/60)%60,'m'
	#if tarTi > 50:
	#	break
		
#print diffLists[0]
#print [x[1] for x in distances[0]]
#print distances[01][03][0],sum(diffLists[01][03])
cPickle.dump(distances,open(fileDrop,'wb'))
#cPickle.dump(diffLists,open(fileDrop+'.sub','wb'))


