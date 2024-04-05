#Lab 2 --Towsifa Akhter
#run with ```python3.12 Lab2-pt1-chi2-2d.py "p0" "p1"``` command because ```root-config --python-version``` gives ```3.12.2```

import ROOT, tdrstyle, sys
import numpy as np

f = ROOT.TFile("Lab2_assignment_v2.root")
event = f.Get("fitdata")

parameter1 = "{}".format(sys.argv[1]) #"p0", "p1", "p2", "p3"
parameter2 = "{}".format(sys.argv[2]) #"p0", "p1", "p2", "p3"

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

xbins = 40
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
p_1 = 1.150
p_2 = 16.530
p_3 = 7.102

#Parameter Test
#'''
v_0 = 1.91
v_1 = 1.21
v_2 = 17.41
v_3 = 8.01
f1 = ROOT.TF1("f1", "[0] + [1]*(x)*exp(-(((x-[2])/[3])**2))", 0, 40)
f1.FixParameter(0, p_0)
f1.FixParameter(1, p_1)
f1.FixParameter(2, p_2)
f1.FixParameter(3, p_3)
f1.SetLineColor(ROOT.kRed)
f1.SetLineStyle(9)
f1.SetMarkerSize(0)
h.Fit("f1")
f1.Draw("same")

#print("what", np.exp(-(((20-p_2)/p_3)**2)))

f2 = ROOT.TF1("f2", "[0] + [1]*(x)*exp(-(((x-[2])/[3])**2))", 0, 40)
f2.FixParameter(0, v_0)
f2.FixParameter(1, v_1)
f2.FixParameter(2, v_2)
f2.FixParameter(3, v_3)
f2.SetLineColor(ROOT.kBlue)
f2.SetLineStyle(8)
f2.SetMarkerSize(0)
h.Fit("f2")
f2.Draw("same")

f3 = ROOT.TF1("f3", "[0] + [1]*(x)*exp(-(((x-[2])/[3])**2))", 0, 40)
f3.FixParameter(0, v_0)
f3.FixParameter(1, 1.11)
f3.FixParameter(2, 17.01)
f3.FixParameter(3, 8.41)
f3.SetLineColor(ROOT.kGreen+2)
f3.SetLineStyle(8)
f3.SetMarkerSize(0)
h.Fit("f3")
f3.Draw("same")

legend_h = ROOT.TLegend(0.6, 0.75, 0.95, 0.95)
legend_h.SetHeader("parameters")
legend_h.AddEntry("f1", "p0 {v0}; p1 {v1}; p2 {v2}; p3 {v3}".format(v0=p_0, v1=p_1, v2=p_2, v3=p_3))
legend_h.AddEntry("f2", "p0 {v0}; p1 {v1}; p2 {v2}; p3 {v3}".format(v0=v_0, v1=v_1, v2=v_2, v3=v_3))
legend_h.AddEntry("f3", "p0 {v0}; p1 {v1}; p2 {v2}; p3 {v3}".format(v0=v_0, v1=1.11, v2=17.01, v3=8.41))
legend_h.SetTextSize(0.)
legend_h.SetBorderSize(0)
legend_h.Draw()
#'''
#end of parameter test

def fit_func(x, d, p1, p2, par1, par2): 
	#x=x invariant mass of bin i #p=the floating parameter 
	#d= number of events in bin i #sigma=uncertainty in bin i = sqrt(f)
	
	if par1 == "p2" and par2== "p3":
		#print(p_0, p_1, p1, p2)
		argument = -pow( ( (x-p1) / p2), 2)
		f = pow(p_0,1) + p_1 * (pow(x,1)) * np.exp(argument)
		chi2 = pow((f-d), 2) / (2*f)
		#print("f: ", argument, f, chi2)

	elif par1 == "p1" and par2 == "p2":
		#print(p_0, p_3)
		argument = -pow( ( (x-p2) / p_3), 2)
		f = pow(p_0,1) + p1 * (pow(x,1)) * np.exp(argument)
		chi2 = pow((f-d), 2) / (2*f)

	elif par1 == "p1" and par2 == "p3":
		#print(p_0, p_2)
		argument = -pow( ( (x-p_2) / p2), 2)
		f = pow(p_0,1) + p1 * (pow(x,1)) * np.exp(argument)
		chi2 = pow((f-d), 2) / (2*f)

	elif par1 == "p0" and par2 == "p1":
		#print(p_2, p_3)
		argument = -pow( ( (x-p_2) / p_3), 2)
		f = pow(p1,1) + p2 * (pow(x,1)) * np.exp(argument)
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
canvas2.SetLogz()

if parameter1=="p2" and parameter2=="p3": 
	x_string = "p_{2}"
	y_string = "p_{3}"
if parameter1=="p1" and parameter2=="p2": 
	x_string = "p_{1}"
	y_string = "p_{2}"
if parameter1=="p1" and parameter2=="p3": 
	x_string = "p_{1}"
	y_string = "p_{3}"
if parameter1=="p0" and parameter2=="p1": 
	x_string = "p_{0}^{2}"
	y_string = "p_{1}"

chi2_hist = ROOT.TH2F("chi2_hist", "chi2_hist", 200, 0.01, 20, 200, 0.01, 20) 

chi2_xAxis = chi2_hist.GetXaxis()
chi2_xAxis.SetTitleOffset(1)
chi2_xAxis.SetTitleSize(0.05)
chi2_xAxis.SetTitle("Floating parameter, {p}".format(p=x_string))

chi2_yAxis = chi2_hist.GetYaxis()
chi2_yAxis.SetTitleOffset(1.2)
chi2_yAxis.SetTitleSize(0.05)
chi2_yAxis.SetTitle("Floating parameter, {p}".format(p=y_string))

chi2_zAxis = chi2_hist.GetZaxis()
chi2_zAxis.SetTitleOffset(1.5)
chi2_zAxis.SetTitleSize(0.04)
chi2_zAxis.SetTitle("#chi^{2}_{i,j}")


p1_list = []
p2_list = []
chi2_list = []

#n = 0
for p1 in np.arange(0.01, 20.01, 0.1):
	#n+=1
	for p2 in np.arange(0.01, 20.01, 0.1):
		chi2_sum = 0
		for j in range(xbins):

			#print(j, h.GetBinContent(j))
			if h.GetBinContent(j) != 0:
				chi2_sum += fit_func(j, h.GetBinContent(j), p1, p2, parameter1, parameter2)
			#print("I'm not a smartass, chi2 = ", chi2_sum)

		chi2_hist.Fill(p1, p2, chi2_sum)
		p1_list.append(p1)
		p2_list.append(p2)
		chi2_list.append(chi2_sum)


#chi2_xAxis.SetRangeUser(0.01,3)
#chi2_yAxis.SetRangeUser(6,12)
chi2_hist.SetMarkerSize(0)
chi2_hist.Draw("colz")

min_chi2 = min(chi2_list)
min_p1 = p1_list[chi2_list.index(min_chi2)]
min_p2 = p2_list[chi2_list.index(min_chi2)]

#'''
#The following few lines are to calculate theuncertainty in the parameter. allowing χ2 to vary by ±1 from
#the minimum corresponds to the 68% C.L. on the parameter measurements. <---Following this idea, I am looking
#for chi2 values that are +1 or -1 froom the minimum chi2 value. But simply adding/subtracting 1 does not mean
#the value actually exists. So I am specifying to look for values that match 0.1% in the list to what I want.
#and then finding the corresponding parameter values to find the uncertainty in the parameter.

min_chi2_list_left = []
min_chi2_list_right = []
for idx, k in enumerate(chi2_list):
	if idx < chi2_list.index(min_chi2):
		if (min_chi2+1)-k < 0.01:
			min_chi2_list_left.append(k)
	if idx > chi2_list.index(min_chi2):
		if (min_chi2+1)-k < 0.01:
			min_chi2_list_right.append(k)

uncertainty_chi_left = min(min_chi2_list_left)
uncertainty_chi_right = min(min_chi2_list_right)
uncertainty_par1_left = p1_list[chi2_list.index(uncertainty_chi_left)]
uncertainty_par1_right = p1_list[chi2_list.index(uncertainty_chi_right)]
uncertainty_par2_left = p2_list[chi2_list.index(uncertainty_chi_left)]
uncertainty_par2_right = p2_list[chi2_list.index(uncertainty_chi_right)]


print("chi2 -1, +1, par1 par2 left, right ", uncertainty_chi_left-min_chi2, uncertainty_chi_right-min_chi2, uncertainty_par1_left-min_p1, uncertainty_par1_right-min_p1, uncertainty_par2_left-min_p2, uncertainty_par2_right-min_p2)

#'''



legend = ROOT.TLegend(0.5, 0.9, 0.9, 1.01)
legend.AddEntry("chi2_hist", "#chi_{min}^{2} "+"{m}".format(m=round(min_chi2, 3) )
				+" at {p1}:{p2}=".format(p1=x_string, p2=y_string)+"{n1}:{n2}".format(n1=round(min_p1, 3), n2=round(min_p2, 3) ) )
legend.SetTextSize(0.)
legend.SetBorderSize(0)
legend.Draw()

line1 = ROOT.TLine(min_p1, 0.01, min_p1, 20)
line1.SetLineColor(ROOT.kRed)
line2 = ROOT.TLine(0.01, min_p2, 20, min_p2)
line2.SetLineColor(ROOT.kRed)
line1.Draw()
line2.Draw()

#'''
latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextAngle(0)
latex.SetTextColor(ROOT.kBlack)
latex.SetTextFont(42)
latex.SetTextSize(0.3*canvas.GetTopMargin())
latex.SetTextAlign(32)
latex.DrawLatex(0+4.0*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.5*canvas.GetTopMargin(), 
	"f_{i} = p_{0}^{2} + p_{1} x_{i}^{2} exp(#frac{-(x_{i}-p_{2})^{2}}{p_{3}^{2}})")


canvas.SaveAs("plots/part1/part1_2D_chi2_testingCalculatedFits.png")
#canvas2.SaveAs("plots/part1/part1_2D_chi2_{p1}{p2}_centered.png".format(p1=parameter1, p2=parameter2))
canvas2.SaveAs("plots/part1/part1_2D_chi2_{p1}{p2}_binwith0removed.png".format(p1=parameter1, p2=parameter2))
