import sys
import pickle
import os.path

# OMIM stuff
from urllib2 import Request, urlopen, URLError
import ast
apiKey='ENTER YOUR KEY HERE'

# LiftOver stuff as OMIM is in hg38 and we're on hg19
from pyliftover import LiftOver
lo = LiftOver('hg19', 'hg38')

dataSets=[]
sampleNames=[]

for sample in sys.argv[2:]:
	thisSample=[]
	sampleNames.append(sample.split('/')[-1])
	for i in range(1,23):
		nextFile=sample+'.'+str(i)+'.pickle'
		if os.path.isfile(nextFile):
			print 'Loading:',nextFile
			thisSample.extend(pickle.load(open(nextFile, 'r')))
		else:
			print 'Missing:',nextFile
	dataSets.append(thisSample)

fileDrop=sys.argv[1]

from xlwt import Workbook, easyxf, Formula
with open(fileDrop+'.xls', 'w') as toXlsFile:
	wb = Workbook()
	ulstyle = easyxf('font: underline single')
	for sampleI,sampleName in enumerate(sampleNames):
		ws = wb.add_sheet(sampleName)
		r = 0
		c = 0
	
		ws.col(len(sampleNames)-1+1).width=256*16
		ws.col(len(sampleNames)-1+4).width=256*26
		ws.col(len(sampleNames)-1+6).width=256*26
	
		xlsCols=sampleNames[:]
		xlsCols.pop(sampleI)
		xlsCols.extend(['SortThing','DupDel','EffectSize','Probes','Unreliable Probes','MinSizeKbp','MinLoc','MaxSizeKbp','MaxLoc','Phenotype','Description'])
		style_string = "font: bold on; borders: bottom dashed"
		style = easyxf(style_string)
		for header in xlsCols:
			ws.write(r,c,header,style=style)
			c+=1
	
		lastCall=[]
		for call in dataSets[sampleI]:
			if call['median'] == 0:
				print call
			#continue
			#print call
			if call == lastCall:
				continue
			lastCall=call

			if call['filteredPost'] != 0:
				continue

			r+=1
			c=1
			
			for sampleJ,otherSampleName in enumerate(sampleNames):
				if sampleJ != sampleI:
					maxOverlap=0
					for otherCall in dataSets[sampleJ]:
						if call['chromosome'] == otherCall['chromosome'] and \
							call['minEnd'] > otherCall['minStart'] and \
							call['minStart'] < otherCall['minEnd']:
								lastStart=max(call['minStart'],otherCall['minStart'])
								firstEnd=min(call['minEnd'],otherCall['minEnd'])
								overLap=firstEnd-lastStart
								myLength=call['minEnd']-call['minStart']
								relativeOverlap=float(overLap)/myLength
								maxOverlap=max(maxOverlap,relativeOverlap)
					ws.write(r,c,maxOverlap)
					c+=1
			
			upDown='Dup'
			if call['mean'] < 0:
				upDown = 'Del'
		
			ws.write(r,c,upDown)
			c+=1
		
			ws.write(r,c,call['mean'])
			c+=1
					
			ws.write(r,c,call['probes'])
			c+=1
			
			ws.write(r,c,call['nonOcc'])
			c+=1
		
			ws.write(r,c,(call['minEnd']-call['minStart'])/1000)
			c+=1
		
			ucscLink = 'HYPERLINK("https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr%s'+':'+'%s'+'-'+'%s";"%s")'
			ucscLink = ucscLink % (call['chromosome'],call['minStart'],call['minEnd'],'chr'+str(call['chromosome'])+':'+str(call['minStart'])+'-'+str(call['minEnd']))
			ws.write(r,c,Formula(ucscLink),ulstyle)
			c+=1
		

			ws.write(r,c,(call['maxEnd']-call['maxStart'])/1000)
			c+=1
		
			ucscLink = 'HYPERLINK("https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr%s'+':'+'%s'+'-'+'%s";"%s")'
			ucscLink = ucscLink % (call['chromosome'],call['maxStart'],call['maxEnd'],'chr'+str(call['chromosome'])+':'+str(call['maxStart'])+'-'+str(call['maxEnd']))
			ws.write(r,c,Formula(ucscLink),ulstyle)
			c+=1
		
			
			worstCase=0
			geneList=[]
			for gene in call['genes']:
				print '\n\n'+gene
				url='http://api.europe.omim.org/api/entry/search?search=gene_symbol:'+gene+'+AND+prefix:*&include=geneMap&format=python&apiKey='+apiKey
				#http://api.europe.omim.org/api/entry/search?search
				request = Request(url)
				#print url

				try:
					response = urlopen(request)
					kittens = response.read()
					omimObject = ast.literal_eval(kittens)
					phenoType=''
					subCase=0
					mimNumber=0
					for stepOne in omimObject['omim']['searchResponse']['entryList']:
						#gene=''
						#print omimObject['omim']['searchResponse']['entryList']
						subCase=0
						mimNumber=0
						if 'mimNumber' in stepOne['entry']:
							mimNumber=stepOne['entry']['mimNumber']
							print mimNumber
						if 'geneMap' in stepOne['entry']:
							if 'phenotypeMapList' in stepOne['entry']['geneMap']:
								#gene=str(stepOne['entry']['geneMap']['geneSymbols'])
								for stepTwo in stepOne['entry']['geneMap']['phenotypeMapList']:
									worstCase=max(worstCase,stepTwo['phenotypeMap']['phenotypeMappingKey'])
									subCase=max(subCase,stepTwo['phenotypeMap']['phenotypeMappingKey'])
									#gene=gene+', '+stepTwo['phenotypeMap']['phenotype']
									phenoType=phenoType+'; '+stepTwo['phenotypeMap']['phenotype']
						#if gene!='':
						#print gene
					geneList.append([subCase,mimNumber,gene,phenoType])
					print worstCase

				except URLError, e:
					print 'Something is wrong with OMIM! Got an error code:', e
			print geneList
			ws.write(r,c,worstCase)
			c+=1
			
			#geneList.sort()
			#print sorted(geneList, reverse=True)
			#continue
			for gene in sorted(geneList, reverse=True):
				#print gene,str(gene[1]),str(gene[2])
				if gene[1] != 0:
					omimLink = 'HYPERLINK("http://europe.omim.org/entry/%i";"%s")'
					datStringYo=str(gene[0])+": "+gene[2]+""+gene[3]
					omimLink = omimLink % (gene[1],datStringYo[:255])
					ws.write(r,c,Formula(omimLink),ulstyle)
				else:
					#ws.write(r,c,'?: '+gene[2])
					#gene_symbol:UBR5+AND+prefix:*
					omimLink = 'HYPERLINK("http://europe.omim.org/search?search=gene_symbol:%s+AND+prefix:*";"%s")'
					omimLink = omimLink % (gene[2],'?: '+gene[2])#'chr'+str(tChrom)+':'+str(posStart)+'-'+str(posEnd))
					ws.write(r,c,Formula(omimLink),ulstyle)
				#ws.write(r,c,gene)
				c+=1
			
			sortScore=(call['probes']-call['nonOcc'])/float(5-worstCase)
			if sortScore<0:
				sortScore=(call['probes']-call['nonOcc'])/float(1+worstCase)
			ws.write(r,0,sortScore)

	wb.save(toXlsFile)
