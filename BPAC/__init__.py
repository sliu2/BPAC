from . import _version
__version__ = _version.__version__

import os
import numpy as np
import pysam
import pyDNase
from pybedtools import BedTool as BT
import pyBigWig
import pandas as pd

def getBigWigMean(regions, bigwig_file, non_nan):
    bw=pyBigWig.open(bigwig_file)
    Profile=[]
    if non_nan ==1: #average over non_nan region
        for region in regions: 
            tmp=bw.stats(str(region.chrom), region.start,region.stop)  # nan considered missing
            if tmp[0]==None:                                           # average over non_nan region
               tmp[0]=0
            Profile.append(tmp[0])
    else:           #average over whole region, default
        for region in regions: 
            values=bw.values(str(region.chrom),region.start,region.stop) # nan considered as 0
            Profile.append(np.mean(np.nan_to_num(values))) #       #average over the whole region
    return Profile

def getBamReadMean(regions, bam_file, non_nan): # bam file should be samtools indexed
    samfile=pysam.AlignmentFile(bam_file,"rb")
    if os.path.isfile(bam_file+'.reads'):  #total number of reads in a genome
        with open(bam_file+'.reads') as f:
            reads=int(f.readline().rstrip())
    else:                       # if genomic reads not provided, normalize with the given bam file
        print 'Please wait, getting total number of (properly paired) mapped reads...'
        print 'assuming bam file is for whole genome region'
        flags=pysam.flagstat(bam_file)
        for flag in flags:
            if (flag.split(' ')[3]=='properly' and flag.split(' ')[4]=='paired'):
                reads=(int(flag.split(' ')[0])+int(flag.split(' ')[2]))
                break
    reads=reads/1000000.0
    print reads, ' millions reads'
    Profile=[]
    if non_nan ==1: #average over non_nan region
        for region in regions:
            sumRegion=0
            sumNonNan=0
            iter=samfile.pileup(str(region.chrom),region.start, region.stop)
            for x in iter:
                if (x.reference_pos<region.stop) and (x.reference_pos>=region.start):
                    #sumRegion+=x.nsegments
                    #if not x.nsegments ==0:
                    #    sumNonNan+=1
                    npp=0
                    for y in x.pileups:
                        if y.alignment.is_proper_pair: # counting only properly paired
                            sumRegion+=1
                            npp+=1
                    if npp>0:
                        sumNonNan+=1
            Profile.append((sumRegion/reads)/sumNonNan)
    else:           #average over whole region, default
        for region in regions:
            sumRegion=0
            iter=samfile.pileup(str(region.chrom),region.start, region.stop)
            for x in iter:
                if (x.reference_pos<region.stop) and (x.reference_pos>=region.start):
                    #sumRegion+=x.nsegments
                    for y in x.pileups:
                        if y.alignment.is_proper_pair:
                            sumRegion+=1
            Profile.append((sumRegion/reads)/(region.stop-region.start))
    return Profile

def getBamCutMean(regions, bam_file):
    cuts=pyDNase.BAMHandler(bam_file)
    Profile=[]
    for region in regions:
        cut=cuts[str(region.chrom)+","+str(region.start)+","+str(region.stop)+",+"]
        Profile.append(sum(cut['+']+cut['-'])*1.0/(region.stop-region.start))
    return Profile

def genAttribute(regions, read_bw_file, cut_bw_file, bam_file, TSSbedFile,chromsizeFile,phastConsFile,label):
# if read_bw_file/cut_bw_file non empty, use them first
    feature=[]
    idxset=[0,1,2,4] #chromesome, start, end, signal(PWM)
    
    featuretmp=[x.split('\t') for x in str(regions).rstrip().split('\n')]
    pdbed=pd.DataFrame(featuretmp)
    feature=pdbed[pdbed.columns[idxset]]
    feature.columns = ['chrom', 'start','end','PWM']
    
    print 'generating DHS/ATAC phastCons feature'
#    feature.loc[:,'phastCons']=pd.Series(getBigWigMean(regions, phastConsFile,0),index=feature.index)
    feature.insert(4,'phastCons',getBigWigMean(regions, phastConsFile,0))
    
    print 'generating DHS/ATAC TSS feature'
    nearbyTSS=[]
    nearby = regions.closest(TSSbedFile, d=True, t='first')
    for region in nearby:
        nearbyTSS.append(region[-1])
    feature.loc[:,'TSS']=pd.Series(nearbyTSS,index=feature.index)

    print 'generating DHS/ATAC reads feature'
    if not (read_bw_file == "0" or  (read_bw_file =="")):
        feature.loc[:,'readProfile']=pd.Series(getBigWigMean(regions, read_bw_file,0),index=feature.index)
    else:
        print 'Please wait, getting raw reading from bam file...'
        feature.loc[:,'readProfile']=pd.Series(getBamReadMean(regions, bam_file,0),index=feature.index)
    
    upregion=regions.flank(l=1, r=0, g=chromsizeFile, pct=True)
    downregion=regions.flank(l=0, r=1, g=chromsizeFile, pct=True)
    
    print 'generating up/down region DHS/ATAC read counts feature'
    if not (read_bw_file == "0" or  (read_bw_file =="")):
        feature.loc[:,'readProfileUp']=pd.Series(getBigWigMean(upregion, read_bw_file,0),index=feature.index)
    else:
        print 'Please wait, getting upstream raw reading from bam file...'
        feature.loc[:,'readProfileUp']=pd.Series(getBamReadMean(upregion, bam_file,0),index=feature.index)
    
    if not (read_bw_file == "0" or  (read_bw_file =="")):
        feature.loc[:,'readProfileDown']=pd.Series(getBigWigMean(downregion, read_bw_file,0),index=feature.index)
    else:
        print 'Please wait, getting downstream raw reading from bam file...'
        feature.loc[:,'readProfileDown']=pd.Series(getBamReadMean(downregion, bam_file,0),index=feature.index)
    
    if not (cut_bw_file == "0" or  (cut_bw_file =="")):
        feature.loc[:,'cutProfile']=pd.Series(getBigWigMean(regions, cut_bw_file,0),index=feature.index)
    else:
        print 'Please wait, getting cut from bam file...'
        feature.loc[:,'cutProfile']=pd.Series(getBamCutMean(regions, bam_file),index=feature.index)
    
    print 'generating up/down region DHS/ATAC cut counts feature'
    if not (cut_bw_file == "0" or  (cut_bw_file =="")):
        feature.loc[:,'cutProfileUp']=pd.Series(getBigWigMean(upregion, cut_bw_file,0),index=feature.index)
    else:
        print 'Please wait, getting upstream cut from bam file...'
        feature.loc[:,'cutProfileUp']=pd.Series(getBamCutMean(upregion, bam_file),index=feature.index)
    
    if not (cut_bw_file == "0" or  (cut_bw_file =="")):
        feature.loc[:,'cutProfileDown']=pd.Series(getBigWigMean(downregion, cut_bw_file,0),index=feature.index)
    else:
        print 'Please wait, getting downstream cut from bam file...'
        feature.loc[:,'cutProfileDown']=pd.Series(getBamCutMean(downregion, bam_file),index=feature.index)
    
    feature.loc[:,'fp_read']=(feature['readProfileUp']+feature['readProfileDown']+1.0)/(feature['readProfile']+1.0)
    feature.loc[:,'fp_cut']=(feature['cutProfileUp']+feature['cutProfileDown']+1.0)/(feature['cutProfile']+1.0)
    feature.loc[:,'label']=int(label)#*feature.shape[0]
    return feature

