#!/usr/bin/env python
# import ROOT in batch mode
import sys
import os
import json
import subprocess
import copy
import tempfile
import time
from commands import getstatusoutput

oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
from optparse import OptionParser
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

# load FWLite C++ librarie
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.gSystem.Load("libDataFormatsEcalDetId.so");
ROOT.gSystem.Load("libPhiSymEcalCalibDataFormats.so");
ROOT.AutoLibraryLoader.enable()

# Load TFile and TTree
from ROOT import TTree, TFile

# getting arrays
#from array import array
import numpy as n

def resetInterval(interval, index):
    interval["index"] = index
    interval["firstRun"] = 0
    interval["firstLumi"] = 0
    interval["unixTimeStart"] = 0
    interval["lastRun"] = 0
    interval["lastLumi"] = 0
    interval["unixTimeEnd"] = 0
    interval["norm"] = 0
    interval["nHit"] = 0
    interval["nLS"] =0
    interval["unixTimeMean"] = 0
    interval["flag"] = ""

def closeInterval(interval):
    interval["unixTimeMean"]=interval["unixTimeStart"]+float(interval["unixTimeMean"])/float(interval["nHit"])

def startInterval(interval, run, lumi, start):
    interval["firstRun"] = run
    interval["firstLumi"] = lumi
    interval["unixTimeStart"] = start

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("", "--debug", dest="debug", action='store_true')
    parser.add_option("", "--save-isolated-intervals", dest="saveIsolatedIntervals", action='store_true')
    parser.add_option("-n", "--max-hit", dest="maxHit", type="int", default=3000000000)
    parser.add_option("-i", "--ls-info-file", dest="ls_info_file", type="string", default="")
    parser.add_option("-o", "--output", dest="output", type="string", default="readMap.root")
    parser.add_option("-t","--max-time", dest="maxTime", type = "int", default=86400)
    parser.add_option("-l","--lumi-file", dest="lumiFile", type = "string", default="")
    (options, args) = parser.parse_args()

    timeMap = {}
    ###--- Load LS info from dataset json file
    with open(options.ls_info_file) as json_file:
        lsInfo = json.load(json_file)
        for time, info in lsInfo.items():
            timeMap[int(time)] = {
                "run"       : info[0],
                "lumi"      : info[1],
                "totHitsEB" : info[2],
                "norm"      : info[3]
                }

    ###--- Certification json
    if options.lumiFile != "":
        with open(options.lumiFile) as json_file:
            goodLumisMap = json.load(json_file)
            if options.debug:
                print(goodLumisMap)
    
    nMaxHits=options.maxHit    
    maxStopTime=options.maxTime

    interval={}

    full_interval_count=0
    isolated_interval_count=0

    currentInterval={}
    resetInterval( currentInterval , 0 )    
    
    # splitting logic
    print("### Start splitting logic")
    for key in sorted(timeMap):
        if options.lumiFile != "" and (timeMap[key]["run"] not in goodLumisMap.keys() or timeMap[key]["lumi"] not in goodLumisMap[timeMap[key]["run"]].keys()):
            print "Skipping lumi: %d:%d not in %s" % (timeMap[key]["run"], timeMap[key]["lumi"], options.lumiFile)
            continue
        
        if currentInterval["nLS"]==0 and currentInterval["unixTimeStart"]==0:
           #start a new interval
            startInterval( currentInterval, timeMap[key]["run"], timeMap[key]["lumi"], key)

        if key-currentInterval["unixTimeStart"]>=maxStopTime and currentInterval["unixTimeStart"] != 0:

            if currentInterval["nHit"] >= nMaxHits/2.:
                # Enough statistics. Closing previous interval 
                if options.debug:
                    print "Closing interval by time condition"

                closeInterval( currentInterval )
                currentInterval["flag"]="S"
                interval[ currentInterval["unixTimeStart" ] ]=dict(currentInterval)
                full_interval_count+=1

            else:
                lastInterval=-1
                if len(interval.keys())>0:
                    lastInterval=sorted(interval.keys())[-1]
                if lastInterval>0:
                    if (currentInterval["unixTimeEnd"]-interval[lastInterval]["unixTimeStart"]<=maxStopTime):
                    #merging with last interval
                        if options.debug:
                            print "Merging interval"
                        closeInterval( currentInterval )
                        interval[lastInterval]["lastRun"]=currentInterval["lastRun"]
                        interval[lastInterval]["lastLumi"]=currentInterval["lastLumi"]
                        interval[lastInterval]["unixTimeEnd"]=currentInterval["unixTimeEnd"]
                        interval[lastInterval]["unixTimeMean"]=(interval[lastInterval]["unixTimeMean"]*interval[lastInterval]["nHit"]+currentInterval["unixTimeMean"]*currentInterval["nHit"])/(float(interval[lastInterval]["nHit"]+currentInterval["nHit"]))
                        interval[lastInterval]["nHit"]+=currentInterval["nHit"]
                        interval[lastInterval]["norm"]+=currentInterval["norm"]
                        interval[lastInterval]["flag"]="M"
                    else:
                        if options.saveIsolatedIntervals:
                            if options.debug:
                                print "Save short interval"
                            closeInterval( currentInterval )
                            currentInterval["flag"]="I"
                            interval[ currentInterval["unixTimeStart" ] ]=dict(currentInterval)
                            full_interval_count+=1
                        else:
                            if options.debug:
                                print "Dropping interval"
                            #dropping interval
                            isolated_interval_count+=1
                else:
                    if options.debug:
                        print "First interval is a short one!"
                    if options.saveIsolatedIntervals:
                        if options.debug:
                            print "Save short interval"
                        closeInterval( currentInterval )
                        currentInterval["flag"]="I"
                        interval[ currentInterval["unixTimeStart" ] ]=dict(currentInterval)
                        full_interval_count+=1
                    else:
                        if options.debug:
                            print "Dropping interval"
                        #dropping interval
                        isolated_interval_count+=1

            # Start a new interval
            resetInterval( currentInterval, full_interval_count )
            startInterval( currentInterval, timeMap[key]["run"], timeMap[key]["lumi"], key)

        currentInterval["lastRun"] = timeMap[key]["run"]
        currentInterval["lastLumi"] = timeMap[key]["lumi"]
        currentInterval["unixTimeEnd"] = key+23.1
        currentInterval["nHit"] += timeMap[key]["totHitsEB"]
        currentInterval["norm"] += timeMap[key]["norm"]
        currentInterval["nLS"] +=1
        currentInterval["unixTimeMean"] += float((key-currentInterval["unixTimeStart"]+11.55)*timeMap[key]["totHitsEB"])

        if currentInterval["nHit"] >= nMaxHits:
            # adding as new interval
            closeInterval( currentInterval )
            currentInterval["flag"]="F"
            interval[ currentInterval["unixTimeStart"] ]=dict(currentInterval)
            full_interval_count+=1
            # resetting for next interval
            resetInterval( currentInterval, full_interval_count )

    interval_number=n.zeros(1,dtype=int)
    hit=n.zeros(1,dtype=long)
    norm=n.zeros(1,dtype=float)
    flag=bytearray(2)
    nLSBranch=n.zeros(1,dtype=int)
    firstRunBranch=n.zeros(1,dtype=int)
    lastRunBranch=n.zeros(1,dtype=int)
    firstLumiBranch=n.zeros(1,dtype=int)
    lastLumiBranch=n.zeros(1,dtype=int)
    unixTimeStartBranch=n.zeros(1,dtype=float)
    unixTimeEndBranch=n.zeros(1,dtype=float)
    unixTimeMeanBranch=n.zeros(1,dtype=float)

    outFile = ROOT.TFile(options.output, "RECREATE")
    if not outFile:
        print "Cannot open outputFile "+options.output

    tree = ROOT.TTree('iov_map', 'IOV map')
    tree.Branch('index', interval_number, 'index/I')
    tree.Branch('flag', flag, 'flag/C')
    tree.Branch('nHit', hit, 'nHit/L')
    tree.Branch('norm', norm, 'norm/D')
    tree.Branch('nLS', nLSBranch, 'nLS/I')
    tree.Branch('firstRun', firstRunBranch, 'firstRun/I')
    tree.Branch('lastRun', lastRunBranch, 'lastRun/I')
    tree.Branch('firstLumi', firstLumiBranch, 'firstLumi/I')
    tree.Branch('lastLumi', lastLumiBranch, 'lastLumi/I')
    tree.Branch('unixTimeStart', unixTimeStartBranch, 'unixTimeStart/D')
    tree.Branch('unixTimeEnd', unixTimeEndBranch, 'unixTimeEnd/D')
    tree.Branch('unixTimeMean', unixTimeMeanBranch, 'unixTimeMean/D')

    for key in sorted(interval):
        interval_number[0]=interval[key]["index"]
        flag[0]=interval[key]["flag"][0]
        hit[0]=interval[key]["nHit"]
        norm[0]=interval[key]["norm"]
        nLSBranch[0]=interval[key]["nLS"]
        firstRunBranch[0]=interval[key]["firstRun"]
        lastRunBranch[0]=interval[key]["lastRun"]
        firstLumiBranch[0]=interval[key]["firstLumi"]
        lastLumiBranch[0]=interval[key]["lastLumi"]
        unixTimeStartBranch[0]=interval[key]["unixTimeStart"]
        unixTimeEndBranch[0]=interval[key]["unixTimeEnd"]
        unixTimeMeanBranch[0]=interval[key]["unixTimeMean"]
        tree.Fill()

    outFile.Write()
    outFile.Close()

    if options.debug:
        for key in sorted(interval):
            print "------------------"
            print "Index: " + str(interval[key]["index"])
            print "nHit: " + str(interval[key]["nHit"])
            print "norm: " + str(interval[key]["norm"])
            print "nLS: " + str(interval[key]["nLS"])
            print "First Run: " + str(interval[key]["firstRun"])
            print "Last Run: " + str(interval[key]["lastRun"])
            print "First Lumi: " + str(interval[key]["firstLumi"])
            print "Last Lumi: " + str(interval[key]["lastLumi"])
            print "Unix Time Start: " + str(interval[key]["unixTimeStart"])
            print "Unix Time End: " + str(interval[key]["unixTimeEnd"])
            print "Unix Time Mean: " + str(interval[key]["unixTimeMean"])
            print "------------------"

    print "====> FULL_INTERVALS:"+str(full_interval_count) + " ISOLATED INTERVALS:" + str(isolated_interval_count)



                  



