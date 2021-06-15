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

    TH1F* hist = new TH1F("a", "N distribution", bin_count, edges);
    TH1F* hist_expectation = new TH1F("b", "2nd dist", bin_count, edges);
    TH1F* hist_correction = new TH1F("c", "3rd dist", bin_count, edges);


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

    hist_correction->GetYaxis()->SetRangeUser(0, 90000);
    

   
    hist_expectation->SetLineColor(kGreen-3);
    hist_correction->SetLineColor(kRed);
    hist->SetLineWidth(4);
    hist_correction->SetLineWidth(3);
    hist_expectation->SetLineWidth(2);
    hist_expectation->SetLineStyle(2);
    //hist_correction->SetLineStyle();

        
    hist_correction->Draw("E");
    hist->Draw("SAME");
    hist_expectation->Draw("SAME");
    //func1->Draw("SAME");

}
