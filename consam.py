import sys
import cPickle
import argparse
import pysam

parser = argparse.ArgumentParser(description='Convert any stream of reads to a pickle file for WISECONDOR, defaults are set for the SAM format',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('bamfile', type=str,
					help='bam file for conversion')
parser.add_argument('bedfile', type=str,
					help='bed file for targets')
parser.add_argument('outfile', type=str,
					help='count file for targets')



args = parser.parse_args()

fileBam		= args.bamfile
fileBed		= args.bedfile
fileOut		= args.outfile
# Prepare the list of chromosomes
chromosomes = dict()
for chromosome in range(1,23):
	chromosomes[str(chromosome)] = []
chromosomes['X'] = []
chromosomes['Y'] = []
chromosomes['M'] = []
targets		= dict()
values		= dict()
curTargets	= []
lastTarget=[0,0,0,0]

with open(fileBed) as inFile:
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
		
samfile = pysam.AlignmentFile(fileBam, "rb")
for chrom in samfile.references:
	print chrom,
	if chrom[3:] in targets:
		curTargets = targets[chrom[3:]]
		print len(curTargets)
		iTarget=iter(curTargets)
		lastTarget=iTarget.next()
		curTarget=iTarget.next()
		lastCount=0
		curCount=0
		index=1
		curVals=[0]*(len(curTargets))
		for read in samfile.fetch(chrom):
			if read.mapping_quality >= 30 and read.is_proper_pair and not read.is_reverse and not read.is_duplicate: # Or the other way around?
				midRead=read.pos+read.isize/2
				while read.pos > (curTarget[0]+curTarget[1])/2:
					curVals[index-1]=lastCount
					index+=1
					lastTarget=curTarget
					#print curTarget
					curTarget=next(iTarget,[1e100,1e100])
					lastCount=curCount
					curCount=0
				#if min(abs(midRead - (curTarget[0]+curTarget[1])/2) , abs(midRead - (lastTarget[0]+lastTarget[1])/2)) <= 250:
				# Pick closest probe
				if abs(midRead - (curTarget[0]+curTarget[1])/2) > abs(midRead - (lastTarget[0]+lastTarget[1])/2):
					# Check if read is fully overlapping probe
					#if lastTarget[0] >= read.pos and lastTarget[1] <= read.pos + read.isize:
					if lastTarget[1] >= read.pos + 20 and lastTarget[0] <= read.pos + read.isize - 20:
						lastCount += 1
						#print midRead,read.pos,read.isize,lastTarget,'L X'
					#elif (midRead - (lastTarget[0]+lastTarget[1])/2) <= 250:
						#print midRead,read.pos,read.isize,lastTarget,'L'
				else:
					# Check if read is fully overlapping probe
					#if curTarget[0] >= read.pos and curTarget[1] <= read.pos + read.isize:
					if curTarget[1] >= read.pos + 20 and curTarget[0] <= read.pos + read.isize - 20:
						curCount += 1
						#print midRead,read.pos,read.isize,curTarget,'C X'
					#elif (midRead - (curTarget[0]+curTarget[1])/2) <= 250:
						#print midRead,read.pos,read.isize,curTarget,'C'

		curVals[index - 1] = lastCount
		if len(curTargets) > index:
			curVals[index] = curCount
		values[chrom] = curVals
cPickle.dump(values, open(fileOut, "wb"))
