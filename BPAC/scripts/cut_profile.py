#!/usr/bin/env python

# Based on 2013 Jason Piper - j.piper@warwick.ac.uk pyDNase script
# add counts on two sides
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse
from clint.textui import progress, puts
import pyDNase

parser = argparse.ArgumentParser(description='Writes WIG file with the cut information based on the regions in reads BED file and the reads in reads BAM file')
parser.add_argument("regions", help="BED file of the regions you want to write wig tracks for")
parser.add_argument("reads", help="The BAM file containing the read data")
parser.add_argument("wig_output", help="Path to write the reads wig track to")
args = parser.parse_args()

reads = pyDNase.BAMHandler(args.reads,caching=True)
regions = pyDNase.GenomicIntervalSet(args.regions)
wigout = open(args.wig_output,"w")

#Required for UCSC upload
print >> wigout, "track type=wiggle_0"

puts("Writing wig tracks...")

for each in progress.bar([item for sublist in sorted(regions.intervals.values()) for item in sorted(sublist, key=lambda peak: peak.startbp)]):
    try:
        prevregionp=str(each.chromosome)+","+str(each.startbp-2)+","+str(each.startbp)+",+"
        prevcuts=reads[prevregionp]
        nextregionp=str(each.chromosome)+","+str(each.endbp+1)+","+str(each.endbp+3)+",+"
        nextcuts=reads[nextregionp]
        pp,pm=prevcuts["+"],prevcuts["-"]
        np,nm=nextcuts["+"],nextcuts["-"]
        cuts=reads[each]
        f,r=cuts["+"], cuts["-"]
        print >> wigout, "fixedStep\tchrom=" + str(each.chromosome) + "\t start="+ str(each.startbp+1) +"\tstep=1"
        counter=0
        for i in f:
            if counter==0: #
                print >> wigout, cuts["+"][counter]+cuts["+"][counter+1]+pp[1]+pp[0]+cuts["-"][counter]+cuts["-"][counter+1]+cuts["-"][counter+2]+pm[1]
            elif (counter==1):
                print >> wigout, cuts["+"][counter]+cuts["+"][counter+1]+cuts["+"][counter-1]+pp[1]+cuts["-"][counter]+cuts["-"][counter+1]+cuts["-"][counter+2]+cuts["-"][0]
            elif (counter==len(f)-2):
                print >> wigout, cuts["+"][counter]+cuts["+"][counter+1]+cuts["+"][counter-1]+cuts["+"][counter-2]+cuts["-"][counter]+cuts["-"][counter+1]+cuts["-"][counter-1]+nm[0]
            elif (counter ==len(f)-1):
                print >> wigout, cuts["+"][counter]+cuts["+"][counter-1]+cuts["+"][counter-2]+np[0]+cuts["-"][counter]+cuts["-"][counter-1]+nm[0]+nm[1]
            else:
                print >> wigout, cuts["+"][counter]+cuts["+"][counter+1]+cuts["+"][counter-1]+cuts["+"][counter-2]+cuts["-"][counter]+cuts["-"][counter+1]+cuts["-"][counter+2]+cuts["-"][counter-1]
            counter=counter+1
    except:
        print "Error: "+str(each.chromosome)+","+str(each.startbp-2)+","+str(each.startbp)+",+, Skipping it"
        continue
