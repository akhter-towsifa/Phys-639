//Lab 1-- Towsifa Akhter

#include <iostream>
#include <string>
#include <math.h>
#include <TH2.h>
#include <TH1.h>
#include <TGraph.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TLegend.h>
#include <TLatex.h>

using namespace std;

//Define functions, parameters, and number of iterations below
#define function(x) 1
//#define function(x) x
//#define function(x) pow(x, 3)
//#define function(x) pow(sin(x), 4)
//#define function(x) 200*exp(-pow((x-3),2)/0.00002) + 2*exp(x-3)


double function_integral_root(double* variable, double* parameter){
	return function(variable[0]);
}

double midPoint_integration(float x, float del_x){
	return function(x)*del_x;
}

double trapezoidal_integration(float x, float del_x){
	float low_x = x-del_x/2; float high_x = x+del_x/2;
	return (function(low_x)+function(high_x))*del_x/2;
}

//double monteCarlo_integration(int n, float lower_limit, float upper_limit){
tuple<float, float, float> monteCarlo_integration(int n, float lower_limit, float upper_limit){
	TRandom random;
	//Initial seed is NOT RANDOM, so results will be the same
	//Setting seed to 0 will give us a new random seed
	random.SetSeed(0);
	double sum = 0.0;
	double sum_square = 0.0;
	TH1F *hist = new TH1F("hist", "hist", n, 0, n);

	for (int i=0; i<n; i++){
		double uniform_dist = random.Uniform(lower_limit, upper_limit);
		//cout << "uniform dist " << uniform_dist << endl;
		hist->Fill(uniform_dist);
		sum += function(uniform_dist);
		sum_square += pow(sum, 2);
	}

	float F_bar_square = pow(sum/n, 2);
	float F_square_bar = sum_square/n;
	float del_F_square = F_square_bar - F_bar_square;
	float uncertainty = sqrt(del_F_square)/sum;

	//return (upper_limit-lower_limit)*sum/n;
	return {(upper_limit-lower_limit)*sum/n, hist->GetRMS(), uncertainty};
}

void methodsComparison()
{
	float lower_limit=0.0; float upper_limit=1; //1 or 4 or M_PI/2
	int arr_iteration[] = {5, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000};
	int bin_size = 10;


	//calculating the analytical integral first using root's integral function TF1::Integral() and corresponding error from TF1::IntegralError
	TF1 analytical_function("function_integral_root", function_integral_root);
	auto integral_analytical_value = analytical_function.Integral(lower_limit, upper_limit); //integrating the function between x=0 and x=1 or M_PI/2 or 4
	std::cout << "integral: " << integral_analytical_value << std::endl;

	FILE *t = fopen("eq0.csv", "w");
	fprintf(t, "iteration, analytical, midpoint, trapezoidal, MC\n");

	//Create a canvas to draw the histogram
    TCanvas *c1 = new TCanvas("c1", "Comparison of Methods", 1200, 800);
    c1->SetTickx(); c1->SetTicky(); c1->SetGrid();

    auto mg = new TMultiGraph();

    auto g0 = new TGraph(); g0->SetMarkerColor(kBlack); g0->SetLineColor(kBlack);
    auto g1 = new TGraph(); g1->SetMarkerColor(kBlue+2); g1->SetLineColor(kBlue+2);
    auto g2 = new TGraph(); g2->SetMarkerColor(kGreen+2); g2->SetLineColor(kGreen+2);
    auto g3 = new TGraph(); g3->SetMarkerColor(kRed); g3->SetLineColor(kRed);
    auto g4 = new TGraph(); g4->SetMarkerColor(kRed); g4->SetLineColor(kRed);
    auto g5 = new TGraph(); g5->SetMarkerColor(kRed); g5->SetLineColor(kRed);


	for (size_t iter=0; iter<sizeof(arr_iteration)/sizeof(int); iter++){
		cout << "iteration: " << arr_iteration[iter] << endl;
		float integration_mid=0.0, integration_trap=0.0, integration_mc=0.0;
		float xi = lower_limit;


		float delta_xi = (upper_limit - lower_limit)/arr_iteration[iter];
		//cout << "delta_xi: " << delta_xi << endl;

		for (int i=1; i <=arr_iteration[iter]; i++){
			xi = i*delta_xi;
			//cout << "xi: " << xi << endl;

			integration_mid +=midPoint_integration(xi, delta_xi);
			integration_trap += trapezoidal_integration(xi, delta_xi);
			//integration_mc = monteCarlo_integration(arr_iteration[iter], lower_limit, upper_limit); 
		}

		auto [mc_integral, rms, uncertainty] = monteCarlo_integration(arr_iteration[iter], lower_limit, upper_limit);

		cout << "integration (mid, trap, mc): " << integration_mid << ", " << integration_trap << ", " << mc_integral << endl;
		fprintf(t, "%i, %f, %f, %f, %f\n", arr_iteration[iter], integral_analytical_value, integration_mid, integration_trap, mc_integral);

		g0->SetPoint(iter, arr_iteration[iter], integral_analytical_value);
		g1->SetPoint(iter, arr_iteration[iter], integration_mid);
		g2->SetPoint(iter, arr_iteration[iter], integration_trap);
		g3->SetPoint(iter, arr_iteration[iter], mc_integral);
		g4->SetPoint(iter, arr_iteration[iter], rms);
		g5->SetPoint(iter, arr_iteration[iter], uncertainty);

	}
	g0->SetTitle("Analytical value");			g0->SetMarkerStyle(33); g0->SetLineStyle(9); g0->SetLineWidth(4);
	g1->SetTitle("MidPoint Integration");		g1->SetMarkerStyle(3);  g1->SetLineStyle(7); g1->SetLineWidth(4);
	g2->SetTitle("Trapezoidal Integration");	g2->SetMarkerStyle(30); g2->SetLineStyle(6); g2->SetLineWidth(4);
	g3->SetTitle("Monte Carlo Integration");	g3->SetMarkerStyle(22); g3->SetLineStyle(2); g3->SetLineWidth(4);

	mg->Add(g0);
	mg->Add(g1);
	mg->Add(g2);
	mg->Add(g3);

	mg->SetTitle("Comparing Different Integration Methods for f(x)=1; Number of Iteration, N; Current Value of the Integral, I(N)");
	//mg->SetTitle("Comparing Different Integration Methods for f(x)=x; Number of Iteration, N; Current Value of the Integral, I(N)");
	//mg->SetTitle("Comparing Different Integration Methods for f(x)=x^{3}; Number of Iteration, N; Current Value of the Integral, I(N)");
	//mg->SetTitle("Comparing Different Integration Methods for f(x)=#frac{x^{4}}{(1+x^{2})^{3}}; Number of Iteration, N; Current Value of the Integral, I(N)");
	//mg->SetTitle("Comparing Different Integration Methods for f(x)=200 e^{(- #frac{(x-3)^{2}}{0.00002})} + 2 e^{(x-3)}; Number of Iteration, N; Current Value of the Integral, I(N)");
	mg->GetYaxis()->SetTitleOffset(1.2); //mg->GetYaxis()->SetTitleSize(0.04);
	mg->GetXaxis()->SetTitleOffset(1);
	mg->GetXaxis()->SetNoExponent();
	mg->GetXaxis()->SetLimits(5, 100000);
	c1->SetLogx();
	mg->Draw("ALP");

	c1->BuildLegend(0.2, 0.7, 0.7, 0.9, "I=<f(x)>_{(0:1)}", "");

	c1->SaveAs("plots/eq0.pdf");

	//RMS plot
	TCanvas *c2 = new TCanvas("c2", "RMS of Monte Carlo Integration", 1200, 800);
    c2->SetTickx(); c2->SetTicky(); c2->SetGrid();
    g4->SetTitle("RMS for Monte Carlo Integration; Number of Iteration, N; RMS Value of the Integral, I_{RMS}(N)");	g4->SetMarkerStyle(27); g4->SetLineStyle(2); g4->SetLineWidth(4);
    g4->GetYaxis()->SetTitleOffset(1.2);
    g4->GetXaxis()->SetTitleOffset(1);
    g4->GetXaxis()->SetNoExponent();
    g4->GetXaxis()->SetLimits(5, 100000);
    c2->SetLogx();
    g4->Draw("AL");

    c2->BuildLegend(0.2, 0.7, 0.7, 0.9, "I_{RMS}=<f(x)>_{(0:1)}", "");
    c2->SaveAs("plots/eq0_RMS.pdf");

	//uncertainty plot
    TCanvas *c3 = new TCanvas("c3", "Uncertainty of Monte Carlo Integration", 1200, 800);
    c3->SetTickx(); c3->SetTicky(); c3->SetGrid();
    g5->SetTitle("Uncertainty for Monte Carlo Integration; Number of Iteration, N; Uncertainty Value of the Integral, I_{unc}(N)");	g5->SetMarkerStyle(27); g5->SetLineStyle(2); g5->SetLineWidth(4);
    g5->GetYaxis()->SetTitleOffset(1.2);
    g5->GetXaxis()->SetTitleOffset(1);
    g5->GetXaxis()->SetNoExponent();
    g5->GetXaxis()->SetLimits(5, 100000);
    c3->SetLogx();
    g5->Draw("AL");

    c3->BuildLegend(0.2, 0.7, 0.7, 0.9, "I_{uncertainty}=<f(x)>_{(0:1)}", "");
    c3->SaveAs("plots/eq0_uncertainty.pdf");

  fclose(t);
}


