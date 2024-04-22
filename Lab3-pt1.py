#Lab 2 --Towsifa Akhter
#run with ```python3.12 Lab3-pt1.py``` command because ```root-config --python-version``` gives ```3.12.2```

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

xbins = 50
xlow = 0
xhigh = 5

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

h0 = ROOT.TH1D("h0", "h0", xbins, xlow, xhigh)
h = ROOT.TH1D("h", "h", xbins, xlow, xhigh)

h0 = f.Get("h1") #f.h1

xAxis = h.GetXaxis()
xAxis.SetTitleOffset(1)
xAxis.SetTitleSize(0.05)
xAxis.SetTitle("Invariant Mass (GeV/c^{2})")

yAxis = h.GetYaxis()
yAxis.SetTitleOffset(0)
yAxis.SetTitleSize(0.05)
yAxis.SetTitle("Events/0.1 GeV/c^{2}")

h = h0.Clone()
h.SetLineWidth(3)
h.SetMarkerSize(0)
h.Draw("HIST")


#Parameter Test
#'''
p_b = 59.91#51.78
p_s = 29.51#16.82
p_alpha = 0.6
p_m0 = 3.0
p_sigma = 0.3

f1 = ROOT.TF1("f1", "[0]*exp(-[1]*(x)) + [2] * exp(-(x-[3])**2/(2*[4]**2 ))", 0, 5)
f1.FixParameter(0, p_b)
f1.FixParameter(1, p_alpha)
f1.FixParameter(2, p_s)
f1.FixParameter(3, p_m0)
f1.FixParameter(4, p_sigma)
f1.SetLineColor(ROOT.kRed)
f1.SetLineStyle(9)
f1.SetMarkerSize(0)
h.Fit("f1")
f1.Draw("same")


legend_h = ROOT.TLegend(0.5, 0.85, 0.95, 0.95)
legend_h.SetHeader("parameters")
legend_h.AddEntry("f1", "B:{b}, #alpha:{a}, S:{s}, m0:{m0}, #sigma:{sigma}".format(b=round(f1.GetParameter(0),3), a=p_alpha, s=round(f1.GetParameter(2),3), m0=p_m0, sigma=p_sigma))
legend_h.SetTextSize(0.)
legend_h.SetBorderSize(0)
legend_h.Draw()
#'''
#end of parameter test

def fit_func(x, d, parameter_b, parameter_s): 

	arg_1 = -0.6*x
	arg_2 = -1*pow((x-3.0), 2)/(2*pow(0.3,2))

	f = parameter_b * np.exp(arg_1) + parameter_s * np.exp(arg_2)

	#x=x invariant mass of bin i 
	#d= number of events in bin i

	ln_mle = d * math.log(f) - f - math.log( math.factorial(round(d) ) )

	return -1*ln_mle


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
#canvas2.SetLogy()

ln_mle_hist = ROOT.TH2F("ln_mle_hist", "ln_mle_hist", 600, 0.01, 60, 600, 0.01, 60) 

ln_mle_xAxis = ln_mle_hist.GetXaxis()
ln_mle_xAxis.SetTitleOffset(1)
ln_mle_xAxis.SetTitleSize(0.05)
ln_mle_xAxis.SetTitle("parameter B")

ln_mle_yAxis = ln_mle_hist.GetYaxis()
ln_mle_yAxis.SetTitleOffset(1.2)
ln_mle_yAxis.SetTitleSize(0.05)
ln_mle_yAxis.SetTitle("parameter S")

ln_mle_zAxis = ln_mle_hist.GetZaxis()
ln_mle_zAxis.SetTitleOffset(1.5)
ln_mle_zAxis.SetTitleSize(0.04)
ln_mle_zAxis.SetTitle("-ln(L)")


b_list = []
s_list = []
ln_mle_list = []

#print("test", h.GetBinContent(35))

for b in np.arange(0.01, 60.01, 0.1):

	for s in np.arange(0.01, 60.01, 0.1):
		ln_mle_sum = 0

		for j in range(xbins+1):

			#print(b, s, j, h.GetBinContent(j))
			#if h.GetBinContent(j) != 0:
			ln_mle_sum += fit_func(j, h.GetBinContent(j), b, s)

		#print(ln_mle_sum)
		ln_mle_hist.Fill(b, s, ln_mle_sum)
		b_list.append(b)
		s_list.append(s)
		ln_mle_list.append(ln_mle_sum)


#ln_mle_xAxis.SetRangeUser(8, 20)
#ln_mle_yAxis.SetRangeUser(0.01, 20)
ln_mle_hist.SetMarkerSize(0)
ln_mle_hist.DrawCopy("colz")
ln_mle_hist.SetContour(3)
ln_mle_hist.Draw("cont3 same")
ln_mle_hist.SetLineColor(ROOT.kRed)

min_ln_mle = min(ln_mle_list)
min_b = b_list[ln_mle_list.index(min_ln_mle)]
min_s = s_list[ln_mle_list.index(min_ln_mle)]

line1 = ROOT.TLine(min_b, 0.01, min_b, 60)
line1.SetLineColor(ROOT.kYellow)
line2 = ROOT.TLine(0.01, min_s, 60, min_s)
line2.SetLineColor(ROOT.kYellow)
line1.Draw()
line2.Draw()

legend = ROOT.TLegend(0.5, 0.9, 0.9, 1.01)
legend.AddEntry("line1", "-ln(L)_{min} "+"{m}".format(m=round(min_ln_mle, 3) )
				+" at B:S = {n1}:{n2}".format(n1=round(min_b, 3), n2=round(min_s, 3) ) )
legend.SetTextSize(0.)
legend.SetBorderSize(0)
legend.Draw()

'''
#The following few lines are to calculate the uncertainty in the parameter for maximum 
#likelihood method. We can use the parabolic approximation (slide 21-23 on https://www.desy.de/~rosem/flc_statistics/data/02_parameter_estimation-A.pdf)
#one solution is uncertainty = sqrt(min_p/n) ;(n= number of events i.e. here, 300)
#parabolic approx.: l(mu-sigma_left) = l(mu+sigma_right) = l_min + 0.5

min_ln_mle_list_left = []
min_ln_mle_list_right = []
for idx, k in enumerate(ln_mle_list):
	if idx < ln_mle_list.index(min_ln_mle):
		if (min_ln_mle+0.5)-k < 0.01:
			min_ln_mle_list_left.append(k)
	if idx > ln_mle_list.index(min_ln_mle):
		if (min_ln_mle+0.5)-k < 0.01:
			min_ln_mle_list_right.append(k)

uncertainty_ln_mle_left = min(min_ln_mle_list_left)
uncertainty_ln_mle_right = min(min_ln_mle_list_right)
uncertainty_b_left = b_list[ln_mle_list.index(uncertainty_ln_mle_left)]
uncertainty_b_right = b_list[ln_mle_list.index(uncertainty_ln_mle_right)]
uncertainty_s_left = s_list[ln_mle_list.index(uncertainty_ln_mle_left)]
uncertainty_s_right = s_list[ln_mle_list.index(uncertainty_ln_mle_right)]


print("ln_mle -1: %f, +1: %f, B left: %f right: %f S left: %f, right: %f" %
	(uncertainty_ln_mle_left-min_ln_mle, uncertainty_ln_mle_right-min_ln_mle, uncertainty_b_left-min_b,
	 uncertainty_b_right-min_b, uncertainty_s_left-min_s, uncertainty_s_right-min_s) 
	)


'''


latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextAngle(0)
latex.SetTextColor(ROOT.kBlack)
latex.SetTextFont(42)
latex.SetTextSize(0.3*canvas.GetTopMargin())
latex.SetTextAlign(32)
latex.DrawLatex(0+4.0*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.5*canvas.GetTopMargin(), 
	"f(m) = B exp(-#alpha m) + S exp(#frac{-(m-m_{0})^{2}}{2 #sigma^{2}})")


canvas.SaveAs("plots/part1/dist.png")
canvas2.SaveAs("plots/part1/part1_2D_mle.png")
#canvas2.SaveAs("plots/part1/part1_2D_mle_{p1}{p2}_centered.png".format(p1=parameter1, p2=parameter2))
