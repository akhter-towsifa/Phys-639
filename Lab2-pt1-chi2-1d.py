#Lab 2 --Towsifa Akhter
#run with ```python3.12 Lab2-pt1-chi2-1d.py "p0"``` command because ```root-config --python-version``` gives ```3.12.2```

import ROOT, tdrstyle, sys
import numpy as np

f = ROOT.TFile("Lab2_assignment_v2.root")
event = f.Get("fitdata")

parameter = "{}".format(sys.argv[1]) #"p0", "p1", "p2", "p3"

ROOT.gROOT.SetBatch(1)
tdrstyle.setTDRStyle()

H_ref = 800
W_ref = 1200
W = W_ref
H = H_ref

T = 0.12*H_ref
B = 0.16*H_ref
L = 0.16*W_ref
R = 0.08*W_ref

xbins = 20
xlow = 0
xhigh = 40

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

h = ROOT.TH1D("h", "h", xbins, xlow, xhigh)

xAxis = h.GetXaxis()
xAxis.SetTitleOffset(1)
xAxis.SetTitleSize(0.05)
xAxis.SetTitle("Invariant Mass (GeV/c^{2})")

yAxis = h.GetYaxis()
yAxis.SetTitleOffset(0)
yAxis.SetTitleSize(0.05)
yAxis.SetTitle("Events/1 GeV/c^{2}")

event.Project("h", "blue")

h.SetLineWidth(3)
h.SetMarkerSize(0)
h.Draw("HIST")

p_0 = 1.201
p_1 = 1.150 #1.150 given and 0.07 guessed
p_2 = 16.530
p_3 = 7.102

#Parameter Test
#'''
v_0 = 4.470
v_1 = 0.39
v_2 = 2.05
v_3 = 39.99
f1 = ROOT.TF1("f1", "[0] + [1]*(x**2)*exp(-(((x-[2])/[3])**2))", 0, 40)
f1.FixParameter(0, p_0)
f1.FixParameter(1, p_1)
f1.FixParameter(2, p_2)
f1.FixParameter(3, p_3)
f1.SetLineColor(ROOT.kRed)
f1.SetLineStyle(9)
f1.SetMarkerSize(0)
h.Fit("f1")
f1.Draw("same")

print("what", np.exp(-(((20-p_2)/p_3)**2)))

f2 = ROOT.TF1("f2", "[0] + [1]*(x**2)*exp(-(((x-[2])/[3])**2))", 0, 40)
f2.FixParameter(0, v_0)
f2.FixParameter(1, v_1)
f2.FixParameter(2, v_2)
f2.FixParameter(3, v_3)
f2.SetLineColor(ROOT.kBlue)
f2.SetLineStyle(8)
f2.SetMarkerSize(0)
h.Fit("f2")
f2.Draw("same")

legend_h = ROOT.TLegend(0.6, 0.75, 0.95, 0.95)
legend_h.SetHeader("parameters")
legend_h.AddEntry("f1", "p0 {v0}; p1 {v1}; p2 {v2}; p3 {v3}".format(v0=p_0, v1=p_1, v2=p_2, v3=p_3))
legend_h.AddEntry("f2", "p0 {v0}; p1 {v1}; p2 {v2}; p3 {v3}".format(v0=v_0, v1=v_1, v2=v_2, v3=v_3))
legend_h.SetTextSize(0.)
legend_h.SetBorderSize(0)
legend_h.Draw()
#'''
#end of parameter test

def fit_func(x, d, p, par): 
	#x=x invariant mass of bin i #p=the floating parameter 
	#d= number of events in bin i #sigma=uncertainty in bin i = sqrt(f)
	
	if par == "p3":
		#print(p_0, p_1, p_2)
		argument = -pow( ( (x-p_2) / p), 2)
		f = pow(p_0,1) + p_1 * (pow(x,2)) * np.exp(argument)
		chi2 = pow((f-d), 2) / (2*f)
		#print("f: ", argument, f, chi2)
		#return chi2

	elif par == "p2":
		#print(p_0, p_1, p_3)
		argument = -pow( ( (x-p) / p_3), 2)
		f = pow(p_0,1) + p_1 * (pow(x,2)) * np.exp(argument)
		chi2 = pow((f-d), 2) / (2*f)
		#return chi2

	elif par == "p1":
		#print(p_0, p_2, p_3)
		argument = -pow( ( (x-p_2) / p_3), 2)
		f = pow(p_0,1) + p * (pow(x,2)) * np.exp(argument)
		chi2 = pow((f-d), 2) / (2*f)
		#return chi2

	elif par == "p0":
		#print(p_1, p_2, p_3)
		argument = -pow( ( (x-p_2) / p_3), 2)
		f = pow(p,1) + p_1 * (pow(x,2)) * np.exp(argument)
		chi2 = pow((f-d), 2) / (2*f)

	return chi2


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

if parameter=="p0": x_string = "p_{0}^{2}"
if parameter=="p1": 
	x_string = "p_{1}"
	canvas2.SetLogy()
if parameter=="p2": x_string = "p_{2}"
if parameter=="p3": x_string = "p_{3}"

chi2_hist = ROOT.TGraph() 

chi2_xAxis = chi2_hist.GetXaxis()
chi2_xAxis.SetTitleOffset(1)
chi2_xAxis.SetTitleSize(0.05)
chi2_xAxis.SetTitle("Floating parameter, {p}".format(p=x_string))

chi2_yAxis = chi2_hist.GetYaxis()
chi2_yAxis.SetTitleOffset(0)
chi2_yAxis.SetTitleSize(0.05)
chi2_yAxis.SetTitle("#chi^{2}_{i}")

i_list = []
chi2_list = []
n = 0
for i in np.arange(0.01, 10, 0.01):
	n+=1
	chi2_sum = 0
	for j in range(xbins):

		#print("I'm a big fat dummy, chi2 = ", chi2_sum)
		#print(j, h.GetBinContent(j))
		if h.GetBinContent(j) != 0:
			chi2_sum += fit_func(j, h.GetBinContent(j), i, parameter)
		#print("I'm not a smartass, chi2 = ", chi2_sum)

	chi2_hist.AddPoint(i, chi2_sum)
	i_list.append(i)
	chi2_list.append(chi2_sum)


chi2_hist.SetLineWidth(3)
chi2_hist.SetMarkerSize(0)
chi2_hist.Draw("AC")


min_chi2 = min(chi2_list)
min_p = i_list[chi2_list.index(min_chi2)]

#'''
#The following few lines are to calculate theuncertainty in the parameter. allowing χ2 to vary by ±1 from
#the minimum corresponds to the 68% C.L. on the parameter measurements. <---Following this idea, I am looking
#for chi2 values that are +1 or -1 froom the minimum chi2 value. But simply adding/subtracting 1 does not mean
#the value actually exists. So I am specifying to look for values that match 0.1% in the list to what I want.
#and then finding the corresponding parameter values to find the uncertainty in the parameter.

min_chi2_list_0 = []
min_chi2_list_1 = []
for idx, k in enumerate(chi2_list):
	if idx < chi2_list.index(min_chi2):
		if (min_chi2+1)-k < 0.01:
			min_chi2_list_0.append(k)
	if idx > chi2_list.index(min_chi2):
		if (min_chi2+1)-k < 0.01:
			min_chi2_list_1.append(k)

uncertainty_chi_0 = min(min_chi2_list_0)
uncertainty_chi_1 = min(min_chi2_list_1)
uncertainty_par_0 = i_list[chi2_list.index(uncertainty_chi_0)]
uncertainty_par_1 = i_list[chi2_list.index(uncertainty_chi_1)]

print("chi2 -1, +1, par -1, +1 ", uncertainty_chi_0, uncertainty_chi_1, uncertainty_par_0-min_p, uncertainty_par_1-min_p)
#print("chi2 -1, +1, par -1, +1 ", uncertainty_chi_0, uncertainty_par_0-min_p)
#'''

legend = ROOT.TLegend(0.5, 0.75, 0.9, 0.85)
#legend.SetHeader("Float parameter")
legend.AddEntry("chi2_hist", "Minimum #chi^{2} "+"{m}".format(m=round(min_chi2, 3) )
				+" at {p}= ".format(p=x_string)+"{n}".format(n=round(min_p, 3) ) )
legend.SetTextSize(0.)
legend.SetBorderSize(0)
legend.Draw()


latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextAngle(0)
latex.SetTextColor(ROOT.kBlack)
latex.SetTextFont(42)
latex.SetTextSize(0.3*canvas.GetTopMargin())
latex.SetTextAlign(32)
latex.DrawLatex(0+4.0*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.5*canvas.GetTopMargin(), 
	"f_{i} = p_{0}^{2} + p_{1} x_{i}^{2} exp(#frac{-(x_{i}-p_{2})^{2}}{p_{3}^{2}})")


canvas.SaveAs("plots/test_binsize_20/part1_1D_testingCalculatedFits.png")
canvas2.SaveAs("plots/test_binsize_20/part1_1D_chi2_{p}.png".format(p=parameter))
