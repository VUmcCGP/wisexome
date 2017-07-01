# WISExome
###### Author: [Roy Straver](https://github.com/rstraver)



## REQUIREMENTS
WISExome was developed and tested using Python2.7. Using any other version may cause errors or faulty results.

The list of packages required and tested versions as reported by `pip freeze`:  
matplotlib==2.0.1  
numpy==1.13.0  
pysam==0.10.0  
scipy==0.19.0  
xlwt==1.2.0  



## PREPARING FOR TESTING
To start testing for aberrations, WISExome needs to learn how genomic regions behave when compared to each other. To do so, reference data should be prepared in the same way as the test sample data later on, using the same genomic reference, etc.

There are 2 files you need to download from other sources:  
> captureregions.bed, a file describing where exactly the regions of interest are.  

> ucscrefseq.bed, a file describing the exact start and stop positions of genes and their exons.  



### File conversion
To convert a .bam file to a file for analysis or training, use the `consam.py` conversion script:
```
python consam.py input/samplename.bam captureregions.bed convert/samplename.hits
```

At this point the data is workable but we can correct for length variations in the target regions:
```
python lennormalize.py convert/samplename.hits captureregions.bed leno/samplename.hits
```



### Reference creation
All reference files should be fed into the reference creation scripts. Move the (length normalized) reference samples (part of the `leno/samplename.hits`) into a separate folder (say `refsamples`) and start the reference-finding script for every combination of chromosomes:  
```
for TARGET in `seq 22 -1 1`
do
  for REF in `seq 1 1 22`
  do
    if [ $TARGET != $REF ]
    then
      python prepref.py refsamples/ refdata/$TARGET.$REF.ref chr$TARGET chr$REF
    fi
  done
done
```

The previous step was split per target chromosome and reference chromosome, now the results have to be combined to make a final selection of reference regions for every region on every target chromosome:
```
python takeref.py refdata/ refout/refname
```



### Determining unreliable target regions
At this point we're just going to act like we're about to test a sample, except we feed the script training samples, use an empty file to fake we have occurrence counts per region and use the `-direct` statement to cut it short:
```
touch emptyfile
for CHROM in X `seq 22 -1 1`
do
  for SAMPLE in `ls refsamples/*.hits`
  do
    # Note: This is tricky when your paths are different, 
    # -d/ -f2 may need to be altered to fit your path
    NAME=`echo $SAMPLE | cut -d/ -f2 | cut -d. -f1`
    echo $NAME
    python test.py \
      $SAMPLE \
      refout/refname.$CHROM \
      occout/$NAME.$CHROM \
      $CHROM \
      captureregions.bed \
      ucscrefseq.bed \
      emptyfile \
      -direct \
        >> occout/all.$CHROM
  done
  # Now we can add the received information together
  cat occout/all.$CHROM \
    | grep 'directCallTag' \
    | awk '{count[$1]++}END{for(j in count) print j,count[j]}' \
    | sort -n > occout/occurrence.$CHROM
done
```



## RUNNING TESTS
To test actual test samples we apply nearly the same script as for the occurrence determination, except now we employ the occurrence files we just created, we feed the script test files, and we may add `-plotfile` to enable plotting results. The `-plotfile` argument is completely optional and can be omitted if speed is preferred.
```
for CHROM in X `seq 22 -1 1`
do
  for SAMPLE in `ls testsamples/*.hits`
  do
    # Note: This is tricky when your paths are different, 
    # -d/ -f2 may need to be altered to fit your path
    NAME=`echo $SAMPLE | cut -d/ -f2 | cut -d. -f1`
    echo $NAME
    python test.py \
      $SAMPLE \
      refout/refname.$CHROM \
      testout/$NAME.$CHROM \
      $CHROM \
      captureregions.bed \
      ucscrefseq.bed \
      occout/occurrence.$CHROM \
      -plotfile
  done
done
```


### OMIM annotation export to excel
The data stored after testing (`testout/$NAME.$CHROM` above) can be used to create a nice Excel-sheet, including OMIM information. Several samples can be added in one go, enabling the script to determine any overlap between them. This will result in a sample per sheet, and every sample gets a column to every other sample with information how much (relatively, 0-1) a call overlaps with a call in another sample. This is useful when trios are sequenced.
> Due to licensing details I currently do not provide an API key for OMIM. You can request one here: https://omim.org/api. Fill your key in between the tickmarks in this line: `apiKey='ENTER YOUR KEY HERE'` near the top of the `excel.py` file.
```
python excel.py outfile sampleA.npz sampleB.npz sampleC.npz
```

## TIPS

Do not use reference data from one lab to test samples from another. Every reference file, laboratory and sequencing machine has its own effect on how read depth per bin behaves. Any results obtained by combining files from different origins are unreliable.