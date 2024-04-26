#Lab 2 --Towsifa Akhter
#run with ```python3.12 Lab3-pt2.py``` command because ```root-config --python-version``` gives ```3.12.2```

import ROOT, tdrstyle, sys, math
import numpy as np

f = ROOT.TFile("lab3.root")
#print(*f.GetListOfKeys())

ROOT.gROOT.SetBatch(1)
tdrstyle.setTDRStyle()

H_ref = 800
W_ref = 950
W = W_ref
H = H_ref

T = 0.12*H_ref
B = 0.16*H_ref
L = 0.12*W_ref
R = 0.18*W_ref

xbins = 10
xlow = 0
xhigh = 10

canvas = ROOT.TCanvas("c1", "c1", 100, 100, W, H)
canvas.SetFillColor(0)
canvas.SetBorderMode(0)
canvas.SetFrameFillStyle(0)
canvas.SetFrameBorderMode(0)
canvas.SetLeftMargin( L/W )
canvas.SetRightMargin( R/W )
canvas.SetTopMargin( T/H )
canvas.SetBottomMargin( B/H )
canvas.SetTickx(0)
canvas.SetTicky(0)
canvas.SetGrid()

h2_0 = ROOT.TH1D("h2_0", "h2_0", xbins, xlow, xhigh)
#h2_c = ROOT.TH1D("h2_c", "h2_c", xbins, xlow, xhigh)
h2_0 = f.Get("h2") #f.h1
#h2_c = h2_0.Clone()

xAxis = h2_0.GetXaxis()
xAxis.SetTitleOffset(1)
xAxis.SetTitleSize(0.05)
xAxis.SetTitle("x")

yAxis = h2_0.GetYaxis()
yAxis.SetTitleOffset(0)
yAxis.SetTitleSize(0.05)
yAxis.SetTitle("Events")

h2_0.SetLineWidth(3)
h2_0.SetMarkerSize(0)
h2_0.Draw("HIST")


#Parameter Test
#'''
param = 10

f1 = ROOT.TF1("f1", "[0]", xlow, xhigh)
f1.SetParameter(0, param)
f1.SetLineColor(ROOT.kRed)
f1.SetLineStyle(9)
f1.SetMarkerSize(0)
h2_0.Fit("f1")
f1.Draw("same")


legend_h = ROOT.TLegend(0.5, 0.85, 0.95, 0.95)
legend_h.SetHeader("parameter")
legend_h.AddEntry("f1", "Constant Parameter: {p}#pm{pe}".format(p=round(f1.GetParameter(0),3), pe=round(f1.GetParError(0),3)) )
legend_h.AddEntry("h2_0", "Total number of events: {n}".format(n=int(h2_0.Integral())))
legend_h.SetTextSize(0.)
legend_h.SetBorderSize(0)
legend_h.Draw()
#'''
#end of parameter test

latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextAngle(0)
latex.SetTextColor(ROOT.kBlack)
latex.SetTextFont(42)
latex.SetTextSize(0.3*canvas.GetTopMargin())
latex.SetTextAlign(32)
#latex.DrawLatex(0+4.0*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.5*canvas.GetTopMargin(), 
#	"f(x) = constant")


canvas.SaveAs("plots/part2/part2_dist_h2.png")

canvas3 = ROOT.TCanvas("c3", "c3", 100, 100, W, H)
canvas3.SetFillColor(0)
canvas3.SetBorderMode(0)
canvas3.SetFrameFillStyle(0)
canvas3.SetFrameBorderMode(0)
canvas3.SetLeftMargin( L/W )
canvas3.SetRightMargin( R/W )
canvas3.SetTopMargin( T/H )
canvas3.SetBottomMargin( B/H )
canvas3.SetTickx(0)
canvas3.SetTicky(0)
canvas3.SetGrid()

h3_0 = ROOT.TH1D("h3_0", "h3_0", xbins, xlow, xhigh)
h3_0 = f.Get("h3") #f.h1
#h = h2.Clone()

xAxis = h3_0.GetXaxis()
xAxis.SetTitleOffset(1)
xAxis.SetTitleSize(0.05)
xAxis.SetTitle("x")

yAxis = h3_0.GetYaxis()
yAxis.SetTitleOffset(0)
yAxis.SetTitleSize(0.05)
yAxis.SetTitle("Events")

h3_0.SetLineWidth(3)
h3_0.SetMarkerSize(0)
h3_0.Draw("HIST")


#Parameter Test
#'''
param = 10

f3 = ROOT.TF1("f3", "[0]", xlow, xhigh)
f3.SetParameter(0, param)
f3.SetLineColor(ROOT.kRed)
f3.SetLineStyle(9)
f3.SetMarkerSize(0)
h3_0.Fit("f3")
f3.Draw("same")


legend_h3 = ROOT.TLegend(0.5, 0.85, 0.95, 0.95)
legend_h3.SetHeader("parameter")
legend_h3.AddEntry("f3", "Constant Parameter: {p}#pm{pe}".format(p=round(f3.GetParameter(0),3), pe=round(f3.GetParError(0),3) ))
legend_h3.AddEntry("h3_0", "Total number of events: {n}".format(n=int(h3_0.Integral())))
legend_h3.SetTextSize(0.)
legend_h3.SetBorderSize(0)
legend_h3.Draw()
#'''
#end of parameter test

latex3 = ROOT.TLatex()
latex3.SetNDC()
latex3.SetTextAngle(0)
latex3.SetTextColor(ROOT.kBlack)
latex3.SetTextFont(42)
latex3.SetTextSize(0.3*canvas.GetTopMargin())
latex3.SetTextAlign(32)
#latex3.DrawLatex(0+4.0*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.5*canvas.GetTopMargin(), 
#	"f(x) = constant")


canvas3.SaveAs("plots/part2/part2_dist_h3.png")

