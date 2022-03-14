#!/usr/bin/env python

#This script combines the histos from the multithreading job
import os
import subprocess
import sys

#-5.4  -5   9      11   0.0639   0.      0.
#-5.4  -5   11     13.5 0.0603   0.      0.
#-5.4  -5   13.5   16.5 0.0577   0.      0.
#-5.4  -5   16.5   19.5 0.0557   0.      0.
#-5.4  -5   19.5   22.5 0.0543   0.      0.
#-5.4  -5   22.5   26   0.0531   0.      0.

#5784.0 0.0450 0.0450 6538.0 0.0445 0.0445


def parse(arguments):
    print arguments
    inputFile = open(arguments,"r")
    outputFileName = arguments[0:arguments.find(".txt")]+"_ShearsTable.txt"
    print outputFileName
    outputFile = open(outputFileName,"w")
    
    iLine = 0
    for line in inputFile:
        etaLow = ""
        etaHigh = ""
        if iLine!=0:
            words = line.split()
            etaLow = words[0]
            etaHigh = words[1]
#            print "%s %s" % (etaLow,etaHigh)
#            for iWord in range(2,len(words):
            ptLow = ""
            ptHigh = ""
            Unc = ""
            for iWord in range(3,len(words)):
                if (iWord%3)==0:
                    if iWord+3 >= len(words):
                        print "%s %s %s %s %s %s %s" % (etaLow,etaHigh,words[iWord],"10000.0",
                                                        words[iWord+1],"0.0","0.0")
                        
                        outputFile.write("%s %s %s %s %s %s %s\n" % (etaLow,etaHigh,words[iWord],
                                                                   "10000.0",words[iWord+1],
                                                                   "0.0","0.0"))
                    else:
                        print "%s %s %s %s %s %s %s" % (etaLow,etaHigh,words[iWord],words[iWord+3],
                                                        words[iWord+1],"0.0","0.0")
                        outputFile.write("%s %s %s %s %s %s %s\n" % (etaLow,etaHigh,words[iWord],
                                                                   words[iWord+3],words[iWord+1],
                                                                   "0.0","0.0"))

        iLine+=1



def main(arguments):
    parse("Summer16_23Sep2016BCDV6_DATA_Uncertainty_AK4PFchs.txt")
    parse("Summer16_23Sep2016EFV6_DATA_Uncertainty_AK4PFchs.txt")
    parse("Summer16_23Sep2016GV6_DATA_Uncertainty_AK4PFchs.txt")
    parse("Summer16_23Sep2016HV6_DATA_Uncertainty_AK4PFchs.txt")
    parse("Summer16_23Sep2016V6_MC_Uncertainty_AK4PFchs.txt")



if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
