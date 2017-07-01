import argparse, cPickle, glob, numpy

parser = argparse.ArgumentParser(description='TBA (or CBA perhaps)',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('refdir', type=str,
					help='files containing read hits per target')
parser.add_argument('bedfile', type=str,
					help='bed file for targets')
parser.add_argument('outdir', type=str,
					help='files containing length normalized read hits per target')
					
args = parser.parse_args()


targets		= dict()
curTargets	= []
lastTarget=[0,0,0,0]

with open(args.bedfile) as inFile:
	for line in inFile:
		target=line.split()[:3]
		if target[0][:3] == 'chr':
			target[0]=target[0][3:]
		target[1]=int(target[1])
		target[2]=int(target[2])
		#target[3]=target[3].split(";")[0]
		target.append(0)
		
		if target[0] != lastTarget[0]:
			targets[lastTarget[0]] = curTargets[:]
			#print len(curTargets)
			curTargets = []
		#curTargets.append(target[1]+((target[2]-target[1])/2))
		curTargets.append([target[1], target[2]])
		lastTarget = target

referenceFiles = glob.glob(args.refdir + '/*.hits')
for refFile in referenceFiles:
	print '\tLoading:\t' + refFile, 
	dumpFile = args.outdir+refFile[len(args.refdir):]
	hitFile =  cPickle.load(open(refFile,'r'))	
	for key in hitFile:
		data=hitFile[key]
		curTarget=targets[key[3:]]
		for i,val in enumerate(data):
			data[i] = float(val)/(curTarget[i][1]-curTarget[i][0])

	#sumCount = float(sum([sum(hitFile[x]) for x in hitFile.keys()]))
	#print numpy.median(hitFile['chrX'])/sumCount > 0.00000175

	print '\tDumping:\t' + dumpFile
	cPickle.dump(hitFile, open(dumpFile, "wb"))
