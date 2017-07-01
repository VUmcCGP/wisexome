import cPickle
import argparse

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

args = parser.parse_args()

refData		= args.refdir
fileDrop	= args.dropfile

#sampleHits = cPickle.load(open(fileHits,'r'))

#reference = dict()
tChroms = [] #[str(x) for x in range(1,23)]
tChroms.append('X')

for tChrom in tChroms:
	curRef = dict()
	for rChrom in range(1,23):
		if tChrom != rChrom:
			print tChrom,rChrom
			newData = cPickle.load(open(refData + '/' + str(tChrom) + '.' + str(rChrom) + '.ref','r'))
			#newData = newData[len(newData)/2:]
			curRef[rChrom]=(newData)
	
	outList=[]
	chrLen=len(curRef[curRef.keys()[0]])
	for i in range(chrLen):
		tmpData=[]

		for rChrom in curRef.keys():
			tmpData.extend([x[0],x[1],rChrom] for x in curRef[rChrom][i])
		tmpData.sort()
		outList.append(tmpData[:100])
	#print outList[1][:5]
	#reference[tChrom] = outList
	#print chrLen,len(outList)
	cPickle.dump(outList,open(fileDrop+'.'+str(tChrom),'wb'))
#cPickle.dump(reference,open(fileDrop,'wb'))
