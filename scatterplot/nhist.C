#include <TH1F.h>

void nhist() {
    const int N1 = 24777;
    const int N2 = 85731;
    const int N3 = 9501;

    // Expectation
    const int N1b = 28751;
    const int N2b = 85731;
    const int N3b = 9854;
    
    // Correction factors
    const double c1 = 1.153;
    const double c3 = 1.109;

    const int N1a = N1*c1;
    const int N2a = 85731;
    const int N3a = N3*c3;


    const float zplus = -0.58740105f;
    const float zminus = 0.58740105f;
    const int bin_count = 3;

    float edges[bin_count + 1] = {-1, zplus, zminus, 1};

    TH1F* hist = new TH1F("hist", "N distribution", bin_count, edges);
    TH1F* hist_expectation = new TH1F("hist_expectation", "2nd dist", bin_count, edges);
    TH1F* hist_correction = new TH1F("hist_correction", "cos #theta*", bin_count, edges);


    //TF1* func1 = new TF1("costheta dist","124373.6816529382*((3*0.311*(1 - x)*(1 - x))/8 + (3*0.687*(1 - x*x))/4 + (3*0.0017*(1 + x)*(1 + x))/8)",-1,1);

    for (int i = 0; i < N1; i++) hist->Fill(-0.9f);
    for (int i = 0; i < N2; i++) hist->Fill( 0.0f);
    for (int i = 0; i < N3; i++) hist->Fill( 0.9f);

    for (int i = 0; i < N1b; i++) hist_expectation->Fill(-0.9f);
    for (int i = 0; i < N2b; i++) hist_expectation->Fill( 0.0f);
    for (int i = 0; i < N3b; i++) hist_expectation->Fill( 0.9f);

    for (int i = 0; i < N1a; i++) hist_correction->Fill(-0.9f);
    for (int i = 0; i < N2a; i++) hist_correction->Fill( 0.0f);
    for (int i = 0; i < N3a; i++) hist_correction->Fill( 0.9f);

    gStyle->SetErrorX(0);
    gStyle->SetOptStat(0);

    hist_correction->GetYaxis()->SetRangeUser(0, 90000);
    
    auto legend = new TLegend(0.775,0.7125,0.975,0.85);
    legend->AddEntry(hist_expectation, "SM Expectation","l");
    legend->AddEntry(hist_correction, "Corrected Distribution","ep");
    legend->AddEntry(hist, "ATLAS Data","ep");
    
    hist_correction->GetXaxis()->SetTitle("cos#theta*");
    hist_correction->GetYaxis()->SetTitle("Events / Bin");
    

    
    
   
    hist_expectation->SetLineColor(kGreen+2);
    hist_expectation->SetLineStyle(2);
    hist_expectation->SetLineWidth(2);
    hist_correction->SetMarkerColor(kRed);
    hist_correction->SetLineColor(kRed);
    hist_correction->SetMarkerStyle(8);    
    hist->SetMarkerStyle(8);
    hist->SetMarkerColor(kBlue);
    hist->SetLineColor(kBlue);
    hist_correction->SetBinError(1,401.4806632);
    hist_correction->SetBinError(3,458.8026947);
    hist_expectation->SetBinError(1,228.7668749);
    hist_expectation->SetBinError(2,6.265269348e-7);
    hist_expectation->SetBinError(3,32.74223966);
    
    TCanvas* canvas1  = new TCanvas("canvas1", "cfit", 80, 80, 1200, 700);
        
    hist_correction->Draw("E1");
    hist->Draw("SAME E1");
    hist_expectation->Draw("SAME E1");
    hist_expectation->ShowBackground()->SetLineColor(kGreen);

    
    legend->Draw();
  
    canvas1->SaveAs("../figures/nhist.pdf");

}
