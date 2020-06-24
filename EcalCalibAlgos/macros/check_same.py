import sys
import os
import csv

import argparse 

from ROOT import gROOT, gStyle, gDirectory
from ROOT import TFile, TCanvas, TGraph, TLine, TH1F,TH2F,TLegend

def doHotSpotList(filename,DEBUG=False):
  spot_lst = []

  with open(filename, 'r') as f:
    reader = csv.reader(f)
    spot_lst= list(reader)
    #for line in reader:
    #    spot_lst.append([int(x) for x in line])

  if DEBUG: print "List of the hot spot crystals", spot_lst
  return spot_lst


#parse arguments
parser = parser = argparse.ArgumentParser()
parser.add_argument('--debug',          action='store_true',        default=False,      help='debugging mode')
parser.add_argument("-p", "--basedir",    action="store",      type=str,         default="/eos/user/f/fcetorel/www/PhiSym/eflow/checkHotSpot_EE_2018_Promptbaseline_4dIOV/",       help="base directory")
parser.add_argument("-s", "--subfold",    action="store",      type=str,         default="",       help="subfolder for the specif case")

args = parser.parse_args()

if (args.debug):
  print "Entering Debugging Mode ... "

plotdir = args.basedir + "/"+args.subfold+"/"
DEBUG=args.debug

hs_list=doHotSpotList("list.csv", DEBUG)
match=0
nomatch=0
for xstal in hs_list:
  ix=xstal[0]
  iy=xstal[1]
  iz=xstal[2]

  filename=plotdir+"Flag_ix_"+ix+"_iy_"+iy+"_iz_"+iz+".png"
  if not os.path.exists(filename): 

    nomatch+=1
    if DEBUG: print '>>> This File does not exist --> '+ filename +'... No match ==', nomatch
    
  else:
    match+=1
    os.system('cp '+filename +" "+ plotdir+"Match_ix_"+ix+"_iy_"+iy+"_iz_"+iz+".png")
    print '>>> This has a match in Flag list :' , filename 
    if DEBUG: print '... Match == ', match

print "# of hotspot xstal are " , len(hs_list)
print "Matched ones are: ", match
print "Non Matched ones are: ", nomatch




