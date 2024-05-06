#Lab 2 --Towsifa Akhter
#run with ```python3.12 Lab3-pt2.py``` command because ```root-config --python-version``` gives ```3.12.2```

import ROOT, tdrstyle, sys, math
import numpy as np

f = ROOT.TFile("lab3.root")
#print(*f.GetListOfKeys())

ROOT.gROOT.SetBatch(1)
tdrstyle.setTDRStyle()

def likelihood_calculator(hist, f):
  ln_mle = 0.0
  for i in range(1, int(hist.GetEntries()+1)):
  	d = hist.GetBinContent(i)
  	ln_mle += -1* d * math.log(f) + f + math.log( math.factorial(round(d) ) )
  	#print(i, d, f, ln_mle)
  return ln_mle

def chi2_calculator(hist, f):
  chi2 = 0.0
  for i in range(1, int(hist.GetEntries()+1)):
    d = hist.GetBinContent(i)
    chi2 += pow((f-d), 2) / f
  return chi2

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

data_ln_mle_h2 = likelihood_calculator(h2_0, f1.GetParameter(0))
data_chi2_h2 = chi2_calculator(h2_0, f1.GetParameter(0))
data_fit_h2 = f1.GetParameter(0)

legend_h = ROOT.TLegend(0.5, 0.85, 0.95, 0.95)
#legend_h.SetHeader("parameter")
legend_h.AddEntry("f1", "Constant Parameter: {p}#pm{pe}".format(p=round(f1.GetParameter(0),3), pe=round(f1.GetParError(0),3)) )
legend_h.AddEntry("h2_0", "Total number of events: {n}".format(n=int(h2_0.Integral())))
legend_h.AddEntry("h2_0", "-ln(L): {l}".format(l=round(data_ln_mle_h2 ,3)) )
legend_h.AddEntry("h2_0", "#chi^{2}: "+"{chi2}".format(chi2=round(data_chi2_h2, 3) ))
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

data_ln_mle_h3 = likelihood_calculator(h3_0, f3.GetParameter(0))
data_chi2_h3 = chi2_calculator(h3_0, f3.GetParameter(0))
data_fit_h3 = f3.GetParameter(0)

legend_h3 = ROOT.TLegend(0.5, 0.85, 0.95, 0.95)
#legend_h3.SetHeader("parameter")
legend_h3.AddEntry("f3", "Constant Parameter: {p}#pm{pe}".format(p=round(f3.GetParameter(0),3), pe=round(f3.GetParError(0),3) ))
legend_h3.AddEntry("h3_0", "Total number of events: {n}".format(n=int(h3_0.Integral())))
legend_h3.AddEntry("h3_0", "-ln(L): {l}".format(l=round(data_ln_mle_h3 ,3)) )
legend_h3.AddEntry("h3_0", "#chi^{2}: "+"{chi2}".format(chi2=round(data_chi2_h3, 3) ))
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


###########################
H_ref = 800
W_ref = 950
W = W_ref
H = H_ref

T = 0.12*H_ref
B = 0.16*H_ref
L = 0.16*W_ref
R = 0.12*W_ref


random = ROOT.TRandom()

def pseudo_data(data_ln_mle, data_chi2, data_fit, p_val_counter_mle, p_val_counter_chi2):
  h_tmp = ROOT.TH1D("h_tmp", "h_tmp", xbins, xlow, xhigh)
  random.SetSeed(0)

  for i in range(1, xbins+1):
    p_data = random.Poisson(data_fit)
    h_tmp.SetBinContent(i, p_data)

  pfit = ROOT.TF1("pfit", "[0]", xlow, xhigh)
  pfit.SetParameter(0, data_fit)
  h_tmp.Fit("pfit")

  pseudo_ln_mle = likelihood_calculator(h_tmp, pfit.GetParameter(0))
  pseudo_chi2 = chi2_calculator(h_tmp, pfit.GetParameter(0))
  h_tmp.Reset()

  if pseudo_ln_mle < data_ln_mle:
    p_val_counter_mle += 1

  if pseudo_chi2 < data_chi2:
    p_val_counter_chi2 += 1

  return pseudo_ln_mle, p_val_counter_mle, pseudo_chi2, p_val_counter_chi2


h_p2 = ROOT.TH1D("h_p2", "h_p2", 100, 10, 50) 
h_p3 = ROOT.TH1D("h_p3", "h_p3", 100, 10, 50)
p_val_h2_mle = 0
p_val_h3_mle = 0

h_p2c = ROOT.TH1D("h_p2c", "h_p2c", 100, 0, 100) 
h_p3c = ROOT.TH1D("h_p3c", "h_p3c", 100, 0, 100)
p_val_h2_chi = 0
p_val_h3_chi = 0

iteration = 1000000
for j in range(iteration):
  pseudo_ln_mle_h2, p_val_h2_mle, pseudo_chi2_h2, p_val_h2_chi = pseudo_data(
    data_ln_mle_h2, data_chi2_h2, data_fit_h2, p_val_h2_mle, p_val_h2_chi)
  h_p2.Fill(pseudo_ln_mle_h2)
  h_p2c.Fill(pseudo_chi2_h2)
  pseudo_ln_mle_h3, p_val_h3_mle, pseudo_chi2_h3, p_val_h3_chi =pseudo_data(
    data_ln_mle_h3, data_chi2_h3, data_fit_h3, p_val_h3_mle, p_val_h3_chi)
  h_p3.Fill(pseudo_ln_mle_h3)
  h_p3c.Fill(pseudo_chi2_h3)


canvas2 = ROOT.TCanvas("c2", "c2", 100, 100, W, H)
canvas2.SetFillColor(0)
canvas2.SetBorderMode(0)
canvas2.SetFrameFillStyle(0)
canvas2.SetFrameBorderMode(0)
canvas2.SetLeftMargin( L/W )
canvas2.SetRightMargin( R/W )
canvas2.SetTopMargin( T/H )
canvas2.SetBottomMargin( B/H )
canvas2.SetTickx(0)
canvas2.SetTicky(0)
canvas2.SetGrid()

hp2_x = h_p2.GetXaxis()
hp2_x.SetTitleOffset(1)
hp2_x.SetTitleSize(0.05)
hp2_x.SetTitle("-ln(L)")

hp2_y = h_p2.GetYaxis()
hp2_y.SetTitleOffset(0)
hp2_y.SetTitleSize(0.05)
hp2_y.SetTitle("Entries")
hp2_y.SetMaxDigits(3)

h_p2.SetLineColor(ROOT.kBlue)
h_p2.SetLineWidth(3)
h_p2.SetMarkerSize(0)
h_p2.SetFillColor(ROOT.kBlue-10)
h_p2.Draw("hist")

line1 = ROOT.TLine(data_ln_mle_h2, 0, data_ln_mle_h2, h_p2.GetMaximum() )
line1.SetLineColor(ROOT.kRed)
line1.SetLineWidth(4)
line1.Draw()

latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextAngle(0)
latex.SetTextColor(ROOT.kBlack)
latex.SetTextFont(42)
latex.SetTextSize(0.3*canvas.GetTopMargin())
latex.SetTextAlign(32)
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.3*canvas.GetTopMargin(), "Entries: {entries}".format(entries = int(h_p2.GetEntries())))
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.7*canvas.GetTopMargin(), "Mean: {mean}".format(mean = round(h_p2.GetMean(),3)))
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-1.1*canvas.GetTopMargin(), "Std Dev: {stddev}".format(stddev = round(h_p2.GetStdDev(),3)))
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-1.5*canvas.GetTopMargin(), "# of bins: 100")
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-2.3*canvas.GetTopMargin(), "p-value: {p}".format(p=round(p_val_h2_mle/iteration, 3)))
latex.SetTextColor(ROOT.kRed)
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-1.9*canvas.GetTopMargin(), "-ln(L)_{data}"+": {l}".format(l=round(data_ln_mle_h2, 2) ))



canvas2.SaveAs("plots/part2/part2_dist_h2_pseudo_mle.png")


canvas4 = ROOT.TCanvas("c4", "c4", 100, 100, W, H)
canvas4.SetFillColor(0)
canvas4.SetBorderMode(0)
canvas4.SetFrameFillStyle(0)
canvas4.SetFrameBorderMode(0)
canvas4.SetLeftMargin( L/W )
canvas4.SetRightMargin( R/W )
canvas4.SetTopMargin( T/H )
canvas4.SetBottomMargin( B/H )
canvas4.SetTickx(0)
canvas4.SetTicky(0)
canvas4.SetGrid()

hp3_x = h_p3.GetXaxis()
hp3_x.SetTitleOffset(1)
hp3_x.SetTitleSize(0.05)
hp3_x.SetTitle("-ln(L)")

hp3_y = h_p3.GetYaxis()
hp3_y.SetTitleOffset(0)
hp3_y.SetTitleSize(0.05)
hp3_y.SetTitle("Entries")
hp3_y.SetMaxDigits(3)

h_p3.SetLineColor(ROOT.kRed)
h_p3.SetLineWidth(3)
h_p3.SetMarkerSize(0)
h_p3.SetFillColor(ROOT.kRed-10)
h_p3.Draw("hist")

line2 = ROOT.TLine(data_ln_mle_h3, 0, data_ln_mle_h3, h_p3.GetMaximum() )
line2.SetLineColor(ROOT.kBlue)
line2.SetLineWidth(4)
line2.Draw()


latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextAngle(0)
latex.SetTextColor(ROOT.kBlack)
latex.SetTextFont(42)
latex.SetTextSize(0.3*canvas.GetTopMargin())
latex.SetTextAlign(32)
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.3*canvas.GetTopMargin(), "Entries: {entries}".format(entries = int(h_p3.GetEntries())))
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.7*canvas.GetTopMargin(), "Mean: {mean}".format(mean = round(h_p3.GetMean(),3)))
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-1.1*canvas.GetTopMargin(), "Std Dev: {stddev}".format(stddev = round(h_p3.GetStdDev(),3)))
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-1.5*canvas.GetTopMargin(), "# of bins: 100")
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-2.3*canvas.GetTopMargin(), "p-value: {p}".format(p=round(p_val_h3_mle/iteration, 3)))
latex.SetTextColor(ROOT.kBlue)
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-1.9*canvas.GetTopMargin(), "-ln(L)_{data}"+": {l}".format(l=round(data_ln_mle_h3, 2) ))


canvas4.SaveAs("plots/part2/part2_dist_h3_pseudo_mle.png")

#####chi2

canvas2c = ROOT.TCanvas("c2c", "c2c", 100, 100, W, H)
canvas2c.SetFillColor(0)
canvas2c.SetBorderMode(0)
canvas2c.SetFrameFillStyle(0)
canvas2c.SetFrameBorderMode(0)
canvas2c.SetLeftMargin( L/W )
canvas2c.SetRightMargin( R/W )
canvas2c.SetTopMargin( T/H )
canvas2c.SetBottomMargin( B/H )
canvas2c.SetTickx(0)
canvas2c.SetTicky(0)
canvas2c.SetGrid()

hp2c_x = h_p2c.GetXaxis()
hp2c_x.SetTitleOffset(1)
hp2c_x.SetTitleSize(0.05)
hp2c_x.SetTitle("#chi^{2}")

hp2c_y = h_p2c.GetYaxis()
hp2c_y.SetTitleOffset(0)
hp2c_y.SetTitleSize(0.05)
hp2c_y.SetTitle("Entries")
hp2c_y.SetMaxDigits(3)

h_p2c.SetLineColor(ROOT.kBlue)
h_p2c.SetLineWidth(3)
h_p2c.SetMarkerSize(0)
h_p2c.SetFillColor(ROOT.kBlue-10)
h_p2c.Draw("hist")

line1c = ROOT.TLine(data_chi2_h2, 0, data_chi2_h2, h_p2c.GetMaximum() )
line1c.SetLineColor(ROOT.kRed)
line1c.SetLineWidth(4)
line1c.Draw()

latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextAngle(0)
latex.SetTextColor(ROOT.kBlack)
latex.SetTextFont(42)
latex.SetTextSize(0.3*canvas.GetTopMargin())
latex.SetTextAlign(32)
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.3*canvas.GetTopMargin(), "Entries: {entries}".format(entries = int(h_p2c.GetEntries())))
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.7*canvas.GetTopMargin(), "Mean: {mean}".format(mean = round(h_p2c.GetMean(),3)))
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-1.1*canvas.GetTopMargin(), "Std Dev: {stddev}".format(stddev = round(h_p2c.GetStdDev(),3)))
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-1.5*canvas.GetTopMargin(), "# of bins: 100")
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-2.3*canvas.GetTopMargin(), "p-value: {p}".format(p=round(p_val_h2_chi/iteration, 3)))
latex.SetTextColor(ROOT.kRed)
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-1.9*canvas.GetTopMargin(), "#chi^{2}_{data}"+": {l}".format(l=round(data_chi2_h2, 2) ))



canvas2c.SaveAs("plots/part2/part2_dist_h2_pseudo_chi2.png")


canvas4c = ROOT.TCanvas("c4c", "c4c", 100, 100, W, H)
canvas4c.SetFillColor(0)
canvas4c.SetBorderMode(0)
canvas4c.SetFrameFillStyle(0)
canvas4c.SetFrameBorderMode(0)
canvas4c.SetLeftMargin( L/W )
canvas4c.SetRightMargin( R/W )
canvas4c.SetTopMargin( T/H )
canvas4c.SetBottomMargin( B/H )
canvas4c.SetTickx(0)
canvas4c.SetTicky(0)
canvas4c.SetGrid()

hp3c_x = h_p3c.GetXaxis()
hp3c_x.SetTitleOffset(1)
hp3c_x.SetTitleSize(0.05)
hp3c_x.SetTitle("#chi^{2}")

hp3c_y = h_p3c.GetYaxis()
hp3c_y.SetTitleOffset(0)
hp3c_y.SetTitleSize(0.05)
hp3c_y.SetTitle("Entries")
hp3c_y.SetMaxDigits(3)

h_p3c.SetLineColor(ROOT.kRed)
h_p3c.SetLineWidth(3)
h_p3c.SetMarkerSize(0)
h_p3c.SetFillColor(ROOT.kRed-10)
h_p3c.Draw("hist")

line2c = ROOT.TLine(data_chi2_h3, 0, data_chi2_h3, h_p3c.GetMaximum() )
line2c.SetLineColor(ROOT.kBlue)
line2c.SetLineWidth(4)
line2c.Draw()


latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextAngle(0)
latex.SetTextColor(ROOT.kBlack)
latex.SetTextFont(42)
latex.SetTextSize(0.3*canvas.GetTopMargin())
latex.SetTextAlign(32)
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.3*canvas.GetTopMargin(), "Entries: {entries}".format(entries = int(h_p3c.GetEntries())))
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.7*canvas.GetTopMargin(), "Mean: {mean}".format(mean = round(h_p3c.GetMean(),3)))
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-1.1*canvas.GetTopMargin(), "Std Dev: {stddev}".format(stddev = round(h_p3c.GetStdDev(),3)))
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-1.5*canvas.GetTopMargin(), "# of bins: 100")
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-2.3*canvas.GetTopMargin(), "p-value: {p}".format(p=round(p_val_h3_chi/iteration, 3)))
latex.SetTextColor(ROOT.kBlue)
latex.DrawLatex(0+7*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-1.9*canvas.GetTopMargin(), "#chi^{2}_{data}"+": {l}".format(l=round(data_chi2_h3, 2) ))


canvas4c.SaveAs("plots/part2/part2_dist_h3_pseudo_chi2.png")
