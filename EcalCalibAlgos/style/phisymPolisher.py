from __future__ import division
import sys
import re
import time
import argparse
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT, os
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

parser = argparse.ArgumentParser(description = 'Polish plots with PhiSym style')
parser.add_argument('inputfile' , default="", help='input file')
parser.add_argument('-k', '--key' , default="", help='key of the object to be restyled')
parser.add_argument('-t', '--type' , default="2D", help='type of plot: 1D / 2D')
parser.add_argument('-d', '--drawOpts' , default="", help='graphical option for 1D plot only')
parser.add_argument('--title' , default="", help='plot title, override histo title')
#parser.add_argument('--rel' , action="store_true", default=False, help='dump phisym ICs correction only')

opts = parser.parse_args()

ROOT.gROOT.Macro(os.path.expanduser('rootlogon.C'))

canvas = ROOT.TCanvas("Canvas_1", "Canvas_1", 800, 450);
ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetFrameBorderSize(0);
ROOT.gStyle.SetFrameBorderMode(0);

p1 = ROOT.TPad("logo_pad",  "", 0.0, 0.85, 0.3, 1.0);
p2 = ROOT.TPad("title_pad", "", 0.25, 0.86, 0.9, 1.0);
p3 = ROOT.TPad("plot_pad",  "", 0.0, 0.0, 1.0, 0.9);
p1.SetFillColor(ROOT.kWhite);
p1.SetFillStyle(0);
p1.SetLeftMargin(0);
p3.SetFillColor(ROOT.kWhite);    
p3.SetFillStyle(0);
p3.SetTopMargin(0.05);
p3.SetRightMargin(0.02);
p3.SetBottomMargin(0.13);
p3.SetLeftMargin(0.08);
p3.SetTicks();
if opts.type == "2D":
    p3.SetRightMargin(0.13);
p3.Draw();
p2.Draw();
p1.Draw();

#---Open file
inFile = ROOT.TFile(opts.inputfile, "READ")
if opts.type == "1D":
    plot = inFile.Get(opts.key)
if opts.type == "2D":
    plot = ROOT.TH2F(inFile.Get(opts.key))


tokens_plot = str.split(plot.GetTitle(), ";")
tokens_opts = str.split(opts.title, ";")
for i in range(4):
    if len(tokens_plot) < i+1:
        tokens_plot.append("")
for i in range(4):
    if len(tokens_opts) > i and tokens_opts[i] != "":
        tokens_plot[i] = tokens_opts[i]
    
title_text = tokens_plot[0]
axis_labels = ";#it{"+tokens_plot[1]+"};#it{"+tokens_plot[2]+"};#it{"+tokens_plot[3]+"}"
    
#---Draw PhiSym logo
fontPath = ROOT.gEnv.GetValue("Root.TTFontPath", "")
arialbi = fontPath+"/arialbi.ttf"
img = ROOT.TImage.Open("blanck.jpg")
img.Tile(400, 100);
img.DrawText(100, 10, "PhiSym", 70, "",
              arialbi, ROOT.TImage.kPlain, "phisym_map.xpm")

p1.cd()
img.Draw();

#---Draw title
title = ROOT.TLatex(0.1, 0.64, "#it{"+title_text+"}")
title.SetTextFont(42);
title.SetTextAlign(13);
title.SetTextSize(0.45);

p2.cd()
title.Draw();

#---Draw plot
plot.SetTitle(axis_labels)
plot.GetXaxis().SetLabelFont(42);
plot.GetXaxis().SetLabelSize(0.045);
plot.GetXaxis().SetTitleFont(42);
plot.GetXaxis().SetTitleSize(0.065);
plot.GetXaxis().SetTitleOffset(0.7);
plot.GetYaxis().SetLabelFont(42);
plot.GetYaxis().SetLabelSize(0.045);
plot.GetYaxis().SetTitleFont(42);
plot.GetYaxis().SetTitleSize(0.065);
plot.GetYaxis().SetTitleOffset(0.6);
if opts.type == "2D":
    plot.GetZaxis().SetLabelFont(42);
    plot.GetZaxis().SetLabelSize(0.045);
    plot.GetZaxis().SetTitleFont(42);
    plot.GetZaxis().SetTitleSize(0.065);
    plot.GetZaxis().SetTitleOffset(0.6);

p3.cd()
if opts.type == "1D":
    plot.Draw(opts.drawOpts)
else:
    plot.Draw("COLZ0")

canvas.Print(opts.key+".eps")
canvas.Print(opts.key+".png")
