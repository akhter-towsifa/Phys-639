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
#include <tuple>

using namespace std;

//Define functions, parameters, and number of iterations below

#define function(x) (200*exp(-pow((x-3),2)/0.00002) + 2*exp(x-3))


double function_integral_root(double* variable, double* parameter){
	return function(variable[0]);
}

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
	cout << "sum: " << sum << endl;
	float F_bar_square = pow(sum/n, 2);
	cout << "F_bar_square: " << F_bar_square << endl;
	float F_square_bar = sum_square/pow(n,2);
	cout << "F_square_bar: " << F_square_bar << endl;
	float del_F_square = F_square_bar - F_bar_square;
	cout << "del_F_square: " << del_F_square << endl;
	float uncertainty = sqrt(del_F_square)/sum;
	cout << "uncertainty: " << uncertainty << endl;
	cout << "RMS value: " << hist->GetRMS() << endl;

	//return {(upper_limit-lower_limit)*sum/n, hist->GetRMS(), uncertainty};
    return {(upper_limit-lower_limit)*sum/n, hist->GetStdDev(), uncertainty};
}

tuple<float, float, float> monteCarlo_integration_leftright(int n, float lower_limit, float upper_limit){
	TRandom random;
	//Initial seed is NOT RANDOM, so results will be the same
	//Setting seed to 0 will give us a new random seed
	random.SetSeed(0);
	double sum = 0.0; double sum1 = 0.0;
	double sum_square = 0.0; double sum_square1 = 0.0;
	TH1F *hist = new TH1F("hist", "hist", n, 0, n);
	TH1F *hist1 = new TH1F("hist1", "hist1", n, 0, n);

	for (int i=0; i<n; i++){
		double uniform_dist = random.Uniform(lower_limit, 2.99);
		double uniform_dist1 = random.Uniform(3.01, upper_limit);
		//cout << "uniform dist " << uniform_dist << endl;
		hist->Fill(uniform_dist);
		sum += function(uniform_dist);
		sum_square += pow(sum, 2);
		hist1->Fill(uniform_dist1);
		sum1 += function(uniform_dist1);
		sum_square1 += pow(sum1, 2);
	}
	cout << "sum: " << sum << endl;
	float F_bar_square = pow(sum/n, 2);
	cout << "F_bar_square: " << F_bar_square << endl;
	float F_square_bar = sum_square/pow(n,2);
	cout << "F_square_bar: " << F_square_bar << endl;
	float del_F_square = F_square_bar - F_bar_square;
	cout << "del_F_square: " << del_F_square << endl;
	float uncertainty = sqrt(del_F_square)/sum;
	cout << "uncertainty: " << uncertainty << endl;

	//return (upper_limit-lower_limit)*sum/n;
	return {((2.99-lower_limit)*sum/n + (upper_limit-3.01)*sum1/n), hist->GetRMS(), uncertainty};
}

//double weighted_monteCarlo_integration(int n, float lower_limit, float upper_limit, float mean_g, float sigma_g){
tuple<float, float, float> weighted_monteCarlo_integration(int n, float lower_limit, float upper_limit, float mean_g, float sigma_g){
	TRandom random;
	//Initial seed is NOT RANDOM, so results will be the same
	//Setting seed to 0 will give us a new random seed
	random.SetSeed(0);
	double sum = 0.0; double sumleft = 0.0; double sumright = 0.0;
	double sum_square = 0.0; double sum_squareleft = 0.0; double sum_squareright = 0.0;
	TH1F *hist = new TH1F("hist", "hist", n, 0, n);
	TH1F *hist1 = new TH1F("hist", "hist", n, 0, n);
	TH1F *hist2 = new TH1F("hist", "hist", n, 0, n);

	for (int i=0; i<n; i++){
		double gaus_dist = random.Gaus(mean_g, sigma_g);
		//cout << "gaus dist " << gaus_dist << endl;
		hist->Fill(gaus_dist);
		sum += function(gaus_dist)/gaus_dist;
		sum_square += pow(sum, 2);
	}
	float F_bar_square = pow(sum/n, 2);
	float F_square_bar = sum_square/pow(n,2);
	float del_F_square = F_square_bar - F_bar_square;
	float uncertainty = sqrt(del_F_square)/sum;

	//return (upper_limit-lower_limit)*sum/n;
	//return {(upper_limit-lower_limit)*sum/n, hist->GetRMS(), hist->GetStdDev()};
	return {(3.01-2.99)*sum/n, hist->GetRMS(), uncertainty};
}


void weighted_MC()
{
	float lower_limit=0; float upper_limit=4; //1 or 4 or M_PI/2
	int arr_iteration[] = {5, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000};


	//calculating the analytical integral first using root's integral function TF1::Integral() and corresponding error from TF1::IntegralError
	TF1 analytical_function("function_integral_root", function_integral_root);
	auto integral_analytical_value = analytical_function.Integral(lower_limit, upper_limit); //integrating the function between x=0 and x=1 or M_PI/2 or 4
	std::cout << "integral: " << integral_analytical_value << std::endl;

	FILE *t = fopen("eq4_weighted.csv", "w");
	fprintf(t, "analytical integral value: %f\n", integral_analytical_value);
	fprintf(t, "iteration, MC, MC_weighted\n");

	//Create a canvas to draw the histogram
    TCanvas *c1 = new TCanvas("c1", "Comparison of MC Methods", 1200, 800);
    c1->SetTickx(); c1->SetTicky(); c1->SetGrid();
    //c1->SetLogx();

    auto mg = new TMultiGraph();

    auto g0 = new TGraph(); g0->SetMarkerColor(kBlack); g0->SetLineColor(kBlack);
    auto g3 = new TGraph(); g3->SetMarkerColor(kRed); g3->SetLineColor(kRed);
    auto g4 = new TGraph(); g4->SetMarkerColor(kViolet-2); g4->SetLineColor(kViolet-2);
    auto g5 = new TGraph(); g5->SetMarkerColor(kViolet-2); g5->SetLineColor(kViolet-2);
    auto g6 = new TGraph(); g6->SetMarkerColor(kViolet-2); g6->SetLineColor(kViolet-2);

	for (size_t iter=0; iter<sizeof(arr_iteration)/sizeof(int); iter++){
		cout << "iteration: " << arr_iteration[iter] << endl;
		float integration_mc=0.0, integration_mc_w=0.0;
		float xi = lower_limit;

		float delta_xi = (upper_limit - lower_limit)/arr_iteration[iter];

		auto [mc_integral, rms, uncertainty] = monteCarlo_integration(arr_iteration[iter], lower_limit, upper_limit);
		auto [mc_integral_leftright, rms_rl, uncertainty_rl] = monteCarlo_integration_leftright(arr_iteration[iter], lower_limit, upper_limit);
		auto [mc_w_integral, rms_w, uncertainty_w] = weighted_monteCarlo_integration(arr_iteration[iter], lower_limit, upper_limit, 3, 1);

		cout << "integration (mc, mc_w): " << mc_integral << ", " << mc_integral_leftright+mc_w_integral << endl;
		fprintf(t, "%i & %f & %f\n", arr_iteration[iter], mc_integral, mc_integral_leftright+mc_w_integral);

		g0->SetPoint(iter, arr_iteration[iter], integral_analytical_value);
		g3->SetPoint(iter, arr_iteration[iter], mc_integral);
		g4->SetPoint(iter, arr_iteration[iter], mc_integral_leftright+mc_w_integral);//integration_mc_w);
		g5->SetPoint(iter, arr_iteration[iter], rms_w);
		g6->SetPoint(iter, arr_iteration[iter], uncertainty_w);


	}
	g0->SetTitle("Analytical value");			g0->SetMarkerStyle(33); g0->SetLineStyle(9); g0->SetLineWidth(4);
	g3->SetTitle("Monte Carlo Integration");	g3->SetMarkerStyle(22); g3->SetLineStyle(2); g3->SetLineWidth(4);
	g4->SetTitle("weighted Monte Carlo Integration");	g4->SetMarkerStyle(23); g3->SetLineStyle(8); g3->SetLineWidth(4);
	
	mg->Add(g0);
	mg->Add(g3);
	mg->Add(g4);

	mg->SetTitle("Comparing MC and weighted MC Integration Methods for f(x)=200 e^{(- #frac{(x-3)^{2}}{0.00002})} + 2 e^{(x-3)}; Number of Iteration, N; Current Value of the Integral, I(N)");
	mg->GetYaxis()->SetTitleOffset(1.2);
	mg->GetXaxis()->SetTitleOffset(1);
	mg->GetXaxis()->SetNoExponent();
	mg->GetXaxis()->SetLimits(5, 100000);
	//c1->SetLogx();
	mg->Draw("ALP");

	c1->BuildLegend(0.2, 0.7, 0.7, 0.9, "I=<f(x)>_{(0:4)}", "");

	c1->SaveAs("plots/eq4_weight.pdf");

	//RMS plot
	TCanvas *c2 = new TCanvas("c2", "RMS of the weighted Monte Carlo Integration", 1200, 800);
    c2->SetTickx(); c2->SetTicky(); c2->SetGrid();
    g5->SetTitle("RMS for weighted Monte Carlo Integration; Number of Iteration, N; RMS Value of the Integral, I_{RMS}(N)");	g4->SetMarkerStyle(27); g4->SetLineStyle(2); g4->SetLineWidth(4);
    g5->GetYaxis()->SetTitleOffset(1.2);
    g5->GetXaxis()->SetTitleOffset(1);
    g5->GetXaxis()->SetNoExponent();
    g5->GetXaxis()->SetLimits(5, 100000);
    //c2->SetLogx();
    g5->Draw("ALP");

    c2->BuildLegend(0.2, 0.7, 0.7, 0.9, "I_{RMS}=<f(x)>_{(0:4)}", "");
    c2->SaveAs("plots/eq4_RMS_weight.pdf");

	//uncertainty plot
    TCanvas *c3 = new TCanvas("c3", "Uncertainty of the weighted Monte Carlo Integration", 1200, 800);
    c3->SetTickx(); c3->SetTicky(); c3->SetGrid();
    g6->SetTitle("Uncertainty for weighted MC Integration; Number of Iteration, N; Uncertainty Value of the Integral, I_{unc}(N)");	g5->SetMarkerStyle(27); g5->SetLineStyle(2); g5->SetLineWidth(4);
    g6->GetYaxis()->SetTitleOffset(1.2);
    g6->GetXaxis()->SetTitleOffset(1);
    g6->GetXaxis()->SetNoExponent();
    g6->GetXaxis()->SetLimits(5, 100000);
    //c3->SetLogx();
    g6->Draw("ALP");

    c3->BuildLegend(0.2, 0.7, 0.7, 0.9, "I_{uncertainty}=<f(x)>_{(0:4)}", "");
    c3->SaveAs("plots/eq4_uncertainty_weight.pdf");


	fclose(t);
}


