import sys
import os
import csv

import argparse 

from ROOT import gROOT, gStyle, gDirectory
from ROOT import TFile, TCanvas, TGraph, TLine, TH1F,TH2F,TLegend
from ROOT import EEDetId


def doHotSpotList(filename,DEBUG=False):
  spot_lst = []

  with open(filename, 'r') as f:
    reader = csv.reader(f)
    spot_lst= list(reader)
    #for line in reader:
    #    spot_lst.append([int(x) for x in line])

  if DEBUG: print "List of the hot spot crystals", spot_lst
  return spot_lst

def doRandomListOLD(tree,step=1,DEBUG=False): #good alternative to the ones using the Ecal class, equivalent method
  xstal_lst = []

  for hashId in range(0,14648,step):
    tree.Draw("ix:iy:iring","hashId=="+str(hashId),"goff")
    ix = tree.GetVal(0)
    iy = tree.GetVal(1)
    iring = tree.GetVal(2)
    if iring[0]<0: iz="-1"
    else: iz="1"
    
    xstal = [str(int(ix[0])),str(int(iy[0])), iz]
    if DEBUG: 
      print "This is the HashId ", hashId, " And ix =", ix[0], " , iy =", iy[0], " , iring =" ,iring[0], "iz =", iz
    xstal_lst.append(xstal)
  

  if DEBUG: print "List of the hot spot crystals", xstal_lst
  return xstal_lst

def doRandomList(step,DEBUG=False):
  xstal_lst = []

  for hashId in range(0,14648,step):
    ix = EEDetId(EEDetId.detIdFromDenseIndex(hashId)).ix()
    iy = EEDetId(EEDetId.detIdFromDenseIndex(hashId)).iy()
    zside = EEDetId(EEDetId.detIdFromDenseIndex(hashId)).zside()

    
    xstal = [str(ix),str(iy), str(zside)]
    if DEBUG: 
      print "This is the HashId ", hashId, " And ix =", ix, " , iy =", iy, " , iz =", zside
    xstal_lst.append(xstal)
  

  if DEBUG: print "List of the hot spot crystals", xstal_lst
  return xstal_lst 
def drawPlot( tree, outputdir, option,thr,consecutive, DEBUG=False):
  subfold=""
  if "hotspot" in option:  
    xstal_list= doHotSpotList("../macros/list.csv") #take the list from Alejandro
    tag="HotSpot"
    print "Drawing crystals from hotspot list"
    #xstal_list=[['28', '54', '-1'],['10','11','1']]

  if "rnd" in option:
    xstal_list = doRandomList(tree,100,DEBUG) # list of random xstal: take a xstal every 100
    tag="Random"
    print "Running on random crystals..."
  
  if "flag" in option:
    print "Flagging Crystals with 3 points over ", thr ,"threshold"
    xstal_list = doRandomList(1000,DEBUG)# list of all xstal
    #xstal_list=[['28', '54', '-1'],['28','53','-1'],['29','45','1']]
    tag="Flag"
    subfold="thr_"+str(thr)
    flagged=0 #counter for flagged events

    if consecutive:
      subfold="thr_"+str(thr)+"_consecutive"
      print "Points have to be consecutive"
  try: 
    os.mkdir(outputdir+subfold) 
  except OSError as error: 
    print(error)


  nxstal=0
  for xstal in xstal_list:
    nxstal+=1

    if nxstal%100==0: 
      print ">>> >> > Processing xstal" , nxstal , " out of " ,14648
      sys.stdout.flush()
    skip=False #always False for rnd and HotSpot
    ix=xstal[0]
    iy=xstal[1]
    iz=xstal[2]
    if "-" in iz: iring="<0"
    else: iring=">0"

    gStyle.SetOptFit(0)
    gStyle.SetOptStat(0)
    gROOT.SetBatch(1)

    c=TCanvas()
    c.cd()
    # ic ratio
    haxis=TH1F("","", 1000 ,1.5242720e+09, 1.5408360e+09)
    tree.SetMarkerStyle(24)
    n1 = tree.Draw("ic_ratio_eflow:avg_time","n_events > 8e6 && ix=="+ix+" && iy=="+iy+" && iring"+iring,"goff")
    ic = tree.GetVal(0)

    if "flag" in option:
      skip=False
      cont = 0
      if DEBUG: print "Option Flag --> Looking for at least 3 points over  threshold"
      for i in range(n1):
        if DEBUG: print "#IOV",i," --> abs(iceflow-1) == ", abs(ic[i]-1)
        if abs(ic[i]-1) > thr: 
          cont +=1
          if DEBUG:   print "Found a point over thr --> cont =", cont
        else:
          if (consecutive): 
            cont=0
            if DEBUG:   print "Consecutive requirment: This point isn't over --> cont =", cont
        if cont==3: 
          flagged +=1
          if DEBUG: print "Draw this event, now are " , flagged , " flagged"
          break
      if cont < 3: 
        skip = True
        if DEBUG: print "So I can Skip this event ..."
    if skip : 
      if DEBUG: print "Skipped"
      continue
  
     
      

    h = TGraph(n1,tree.GetVal(1),tree.GetVal(0))
    if DEBUG:

      print "Drawing ic plot: ", "ic_ratio_eflow:avg_time","n_events > 8e6 && ix=="+ix+" && iy=="+iy+" && iring"+iring,""

    h.SetMarkerColor(633)
    h.SetMarkerStyle(24) #open Circle
  
    #lc 

    n2 = tree.Draw("1/lc:avg_time","n_events > 8e6 && ix=="+ix+" && iy=="+iy+" && iring"+iring,"goff")
    hl = TGraph(n2,tree.GetVal(1),tree.GetVal(0))
    if DEBUG: 

      print "Drawing 1/lc plot: ", "1/lc:avg_time>>hl","n_events > 8e6 && ix=="+ix+" && iy=="+iy+" && iring"+iring,""

    hl.SetMarkerColor(1)
    hl.SetMarkerStyle(24) #open Circle


    c.Clear()
    c.SetGrid()
    haxis.SetTitle(";Time (day/month);#it{Response variation};")
    haxis.SetMarkerStyle(0)
    haxis.GetYaxis().SetRangeUser(0.0, 1.2) 
    haxis.GetXaxis().SetTimeDisplay(1) 
    haxis.GetXaxis().SetTimeFormat("%d/%m%F1970-01-01 00:00:00") 
    #haxis.GetXaxis().SetNdivisions(507)
    haxis.Draw("")
    h.Draw("p same")
    hl.Draw("p same")

    line1 = TLine(1.5242720e+09,1-thr,1.5408360e+09,1-thr)
    line1.SetLineWidth(2)
    line1.SetLineStyle(2)
    line1.SetLineColor(2)
    line1.Draw("same")

    line2 = TLine(1.5242720e+09,1+thr,1.5408360e+09,1+thr)
    line2.SetLineWidth(2)
    line2.SetLineStyle(2)
    line2.SetLineColor(2)
    line2.Draw("same")

    leg = TLegend(0.58,0.5,0.89,0.64);
    leg.AddEntry(h,"Response Variation","p")
    leg.AddEntry(hl,"1/ Laser Correction ","p")
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.Draw("same")

    c.Update()

 

    c.SaveAs(outputdir+subfold+"/"+tag+"_ix_"+ix+"_iy_"+iy+"_iz_"+iz+".png")
    c.SaveAs(outputdir+subfold+"/"+tag+"_ix_"+ix+"_iy_"+iy+"_iz_"+iz+".root")
  if "flag" in option: print ">>> >> > These are the # of crystals flagged: ", flagged





#parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--debug',           action='store_true',           default=False,      help='debugging mode')
parser.add_argument("-l", "--label",     action="store",      type=str,                     help="job label: hotspot, rnd or flag")
parser.add_argument("-t", "--thr",     action="store",      type=float,       default=0.05   ,          help="to be used with flag option, choose a threshold es. 0.05 --> 5%")
parser.add_argument('--consecutive',           action='store_true',           default=False,      help='to be used with flag option, consecutive point only')
parser.add_argument("-o", "--outdir",    action="store",      type=str, default="/eos/user/f/fcetorel/www/PhiSym/eflow/checkHotSpot_EE_2018_Promptbaseline_4dIOV/",       help="output directory")
parser.add_argument("-i", "--inputfile", action="store",      type=str,                     help="input file")



args = parser.parse_args()

if (args.debug):
  print "Entering Debugging Mode ... "


outputdir=args.outdir #"/eos/user/f/fcetorel/www/PhiSym/eflow/history_2018_Promptbaseline/"
filename = args.inputfile
if not os.path.exists(filename): 
  print '>>> This File does not exist -->  '+filename+'... check path ... '
else:
  print '>>> Opening File :' , filename
  inFile = TFile.Open ( filename," READ ")
  tree= inFile.Get("ee")
  drawPlot(tree,outputdir, args.label,args.thr,args.consecutive,args.debug)  
  
  
