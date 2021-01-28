import sys
import os
import csv

import argparse 

from ROOT import gROOT, gStyle, gDirectory
from ROOT import TFile, TCanvas, TGraph, TLine, TF1,TH2F,TLegend
from ROOT import EEDetId
from array import array


#Define new branches that are not stored in the eflow ntuples
#perform the fit of the history of the crystal and save the slope and the chi2
def addNewBranches(tree,outdir, savePlot, DEBUG):

    print ">>> Adding new branches..."
    if savePlot: print ">>> ... and saving the plots too ..."
    # create 1 dimensional float arrays as fill variables, in this way the float
    # array serves as a pointer which can be passed to the branch
    slope  = array('d',[0])
    chi2 = array('d',[0])
    iz = array('i',[0])
 
    # create the branches and assign the fill-variables to them as doubles (D)
    bslope = tree.Branch("slope",  slope,  'slope/D')
    bchi2 = tree.Branch("chi2", chi2, 'chi2/D')
    biz = tree.Branch("iz", iz, 'iz/I')
    

    for  i in range( tree.GetEntries()):
        if i%100==0: 
           print ">>> >> > Processing entry" , i , " out of " ,14648
        if DEBUG and i==10 :  break
        tree.GetEntry(i)
        #print "hashID  ", str(tree.hashId)
        hashedId = tree.hashId
        n1 = tree.Draw("ic_ratio_eflow:avg_time","hashId=="+str(tree.hashId),"goff")
        
 
        g=TGraph(n1,tree.GetVal(1),tree.GetVal(0))
        f =  TF1("f","pol1")
        g.SetMarkerColor(633)
        g.SetMarkerStyle(24) #open Circle
        g.Fit(f,"Q")
        
        if savePlot and i%100==0:
            c1 = TCanvas()
	    c1.SetGrid()
	    g.SetTitle(";Time (day/month);#it{Response variation};")
	    g.GetYaxis().SetRangeUser(0.4, 1.2) 
	    g.GetXaxis().SetTimeDisplay(1) 
	    g.GetXaxis().SetTimeFormat("%d/%m%F1970-01-01 00:00:00") 
	    #haxis.GetXaxis().SetNdivisions(507)
            g.Draw("AP")
            c1.SaveAs(outdir+"/xstal_"+str(hashedId)+".png")    
            c1.SaveAs(outdir+"/xstal_"+str(hashedId)+".root")    
        
        slope[0] = f.GetParameter(1)
        chi2[0] = f.GetChisquare()
        iz[0] = EEDetId(EEDetId.detIdFromDenseIndex(hashedId)).zside()
        
        bslope.Fill()
        bchi2.Fill()
        biz.Fill()      

 
          
def  drawOnly(tree,  outdir, channels, DEBUG ):
    for hashedId in channels:
        n1 = tree.Draw("ic_ratio_eflow:avg_time","hashId=="+str(hashedId),"goff")
         
        g=TGraph(n1,tree.GetVal(1),tree.GetVal(0))
        f =  TF1("f","pol1")
        g.SetMarkerColor(633)
        g.SetMarkerStyle(24) #open Circle
        g.Fit(f,"Q")

        c1 = TCanvas()
        c1.SetGrid()
        g.SetTitle(";Time (day/month);#it{Response variation};")
        g.GetYaxis().SetRangeUser(0.0, 1.2) 
        g.GetXaxis().SetTimeDisplay(1) 
        g.GetXaxis().SetTimeFormat("%d/%m%F1970-01-01 00:00:00") 
        g.Draw("AP")
        c1.SaveAs(outdir+"/xstal_"+str(hashedId)+".png")    
        c1.SaveAs(outdir+"/xstal_"+str(hashedId)+".root")    

def  mktree(oldtree, newfile_name, outdir, savePlot, DEBUG ):
    print ">>> Reading the old tree"
    oldtree.SetBranchStatus("*",0)
    # Activate only the branches we need for EE studies
    for  activeBranchName in  {"n_iovs","ix","iy", "hashId","avg_time", "firstRun", "lastRun","ic_ratio_eflow"} :
         oldtree.SetBranchStatus(activeBranchName, 1)
   
                             
    #Create a new file + a clone of old tree header and copy all events
    
    newfile = TFile(newfilename, "recreate")
    newtree = oldtree.CloneTree(0)                                                 
   
    newtree.CopyEntries(oldtree)
                                    
    if DEBUG:
        print "Copied the branches of the old tree"
        newtree.Print()
    
    addNewBranches(newtree, outdir, savePlot, DEBUG)
   
    #change name to a branch -. compatibility
    newtree.GetBranch("hashId").SetTitle("hashedId/I")
    newtree.GetBranch("hashId").SetName("hashedId")
    #newtree.GetBranch("ic_ratio_eflow").SetTitle(ICname+"/F")
    #newtree.GetBranch("ic_ratio_eflow").SetName(ICname)
 
    if DEBUG: newtree.Print()
    newfile.Write()
  

#Produce the tree to be used for EE hotspots studies
#1) Copy only the relevant info stored in the eflow ntuples from step 4
#2) Add some new variable i.e. the slope and chi2 of the fit performed on the history of each crystal
#parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--debug',           action='store_true',           default=False,      help='debugging mode')
parser.add_argument("-o", "--outdir",    action="store",      type=str, default="~/www/PhiSym/eflow/hotspotStudiesEE/UL2018",       help="output directory for plot")
parser.add_argument('--savePlot',           action='store_true',           default=False,      help='save also some xstals plot')
parser.add_argument('-c', '--crystals'   , nargs='+' , default="0" ,  help='Choose the crystal you want to run on (hashed index)')
parser.add_argument('--drawOnly' ,  action='store_true',           default=False,      help='draw only some xstals in the list user specify')
parser.add_argument("-i", "--inputfile", action="store",      type=str,                     help="input file")


gStyle.SetOptFit(1111)
gStyle.SetOptStat(0)
gROOT.SetBatch(1)


args = parser.parse_args()

if (args.debug):
  print ">>> Entering Debugging Mode ... "
DEBUG = args.debug
savePlot = args.savePlot
outdir=args.outdir

outputdir=args.outdir #"/eos/user/f/fcetorel/www/PhiSym/eflow/checkHotSpot_EE_2018_Promptbaseline_4dIOV/"
filename = args.inputfile
if not os.path.exists(filename): 
  print '>>> This File does not exist -->  '+filename+'... check path ... '
else:
  print '>>> Opening File :' , filename
  newfilename = args.inputfile.replace('history', 'hotspotStudies')
  if DEBUG: print 'This is the new file name: ', newfilename
  inFile = TFile.Open ( filename," READ ")
  tree= inFile.Get("ee")
  if args.drawOnly:
    print '>>> Only Draw Mode: drawing only some interesting channels...'
    drawOnly(tree, outdir, args.crystals, DEBUG)    
     
  else: mktree(tree, newfilename,  outdir,  savePlot , DEBUG )



