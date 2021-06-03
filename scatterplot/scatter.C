#include <TFile.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLegend.h>
#include <TGaxis.h>

TLine* make_diagonal_line(TCanvas* canvas){
    canvas->Update();
    double p1 = fmax(canvas->GetUxmin(), canvas->GetUymin());
    double p2 = fmin(canvas->GetUxmax(), canvas->GetUymax());
    
    TLine *line = new TLine(p1, p1, p2, p2);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kRed);
    line->SetLineWidth(3);
    return line;
}

void plot(){
    TFile* delphes_file = new TFile("delphes.root");
    TH1F* truth_0       = delphes_file->Get<TH1F>("CTStarGen0");
    TH1F* truth_R       = delphes_file->Get<TH1F>("CTStarGenR");
    TH1F* truth_L       = delphes_file->Get<TH1F>("CTStarGenL");
    TH1F* costheta_0    = delphes_file->Get<TH1F>("CTStar0");
    TH1F* costheta_R    = delphes_file->Get<TH1F>("CTStarR");
    TH1F* costheta_L    = delphes_file->Get<TH1F>("CTStarL");

    int nbins = truth_0->GetNbinsX();

    double* x_0       = (double*)malloc(sizeof(double)*nbins);
    double* y_0       = (double*)malloc(sizeof(double)*nbins);
    double* x_error_0 = (double*)malloc(sizeof(double)*nbins);
    double* y_error_0 = (double*)malloc(sizeof(double)*nbins);

    double* x_L       = (double*)malloc(sizeof(double)*nbins);
    double* y_L       = (double*)malloc(sizeof(double)*nbins);
    double* x_error_L = (double*)malloc(sizeof(double)*nbins);
    double* y_error_L = (double*)malloc(sizeof(double)*nbins);

    double* x_R       = (double*)malloc(sizeof(double)*nbins);
    double* y_R       = (double*)malloc(sizeof(double)*nbins);
    double* x_error_R = (double*)malloc(sizeof(double)*nbins);
    double* y_error_R = (double*)malloc(sizeof(double)*nbins);
    
    

    for (int i = 0; i < nbins; i++) {
        double x_entries = truth_0   ->GetEntries();
        double y_entries = costheta_0->GetEntries();
        x_0[i]       = truth_0   ->GetBinContent(i+1);
        y_0[i]       = costheta_0->GetBinContent(i+1);
        x_error_0[i] = sqrt(x_0[i])/x_entries;
        y_error_0[i] = sqrt(y_0[i])/y_entries;
        x_0[i] /= x_entries;
        y_0[i] /= y_entries;
    }

    for (int i = 0; i < nbins; i++) {
        double x_entries = truth_L   ->GetEntries();
        double y_entries = costheta_L->GetEntries();
        x_L[i]       = truth_L   ->GetBinContent(i+1);
        y_L[i]       = costheta_L->GetBinContent(i+1);
        x_error_L[i] = sqrt(x_L[i])/x_entries;
        y_error_L[i] = sqrt(y_L[i])/y_entries;
        x_L[i] /= x_entries;
        y_L[i] /= y_entries;
    }
    
    for (int i = 0; i < nbins; i++) {
        double x_entries = truth_R   ->GetEntries();
        double y_entries = costheta_R->GetEntries();
        x_R[i]       = truth_R   ->GetBinContent(i+1);
        y_R[i]       = costheta_R->GetBinContent(i+1);
        x_error_R[i] = sqrt(x_R[i])/x_entries;
        y_error_R[i] = sqrt(y_R[i])/y_entries;
        x_R[i] /= x_entries;
        y_R[i] /= y_entries;
    }

    TCanvas *c0 = new TCanvas("costheta 0", "c0", 80, 80, 1400, 700);
    c0->Range(0, 0, 1, 1);
    c0->SetTitle("Delphes Scatterplot");

    // double max0, maxL, maxR;
    // max0 = maxL = maxR = 0;
    // //double max0, maxL, maxR;
    // for (int i = 0; i < nbins; i++) {
    //     max0 = fmax(x_0[i] + x_error_0[i], max0);
    //     maxL = fmax(x_L[i] + x_error_L[i], maxL);
    //     maxR = fmax(x_R[i] + x_error_R[i], maxR);
    // }
    // printf("max 0 : %g\n", max0);
    // printf("max L : %g\n", maxL);
    // printf("max R : %g\n", maxR);

    // double min0, minL, minR;
    // min0 = minL = minR = 100.0;
    // //double min0, minL, minR;
    // for (int i = 0; i < nbins; i++) {
    //     min0 = fmin(x_0[i] - x_error_0[i], min0);
    //     minL = fmin(x_L[i] - x_error_L[i], minL);
    //     minR = fmin(x_R[i] - x_error_R[i], minR);
    // }
    // printf("min 0 : %g\n", min0);
    // printf("min L : %g\n", minL);
    // printf("min R : %g\n", minR);

    // L is the widest x-interval so use that to create the axis
    TGraph *gL = new TGraphErrors(nbins, x_L, y_L, x_error_L, y_error_L);
    gL->SetTitle("Delphes Scatterplot");
    gL->SetLineColor(kRed);
    gL->SetLineWidth(2);
    gL->Draw("SAME aep");

    TGraph *g0 = new TGraphErrors(nbins, x_0, y_0, x_error_0, y_error_0);
    g0->SetLineColor(kBlue);
    g0->SetLineWidth(2);
    g0->Draw("SAME e");

    TGraph *gR = new TGraphErrors(nbins, x_R, y_R, x_error_R, y_error_R);
    gR->SetLineWidth(2);
    gR->Draw("SAME e");

    gL->GetXaxis()->SetRangeUser(-1, 1.0);
    gL->GetYaxis()->SetRangeUser(0, 0.207376);


    TLine *line0 = make_diagonal_line(c0);
    line0->Draw();

    TLegend* leg = new TLegend(0.1,0.7,0.3,0.9);
    //leg->SetHeader("cos theta* 0");
    leg->AddEntry(gL,    "costheta* L", "lep");
    leg->AddEntry(g0,    "costheta* 0", "lep");
    leg->AddEntry(gR,    "costheta* R", "lep");
    leg->AddEntry(line0, "Expectation","l r");
    leg->Draw();

    c0->SaveAs("scatterplot.png");


/*
    TCanvas *cL = new TCanvas("costheta L", "cL", 80, 80, 700, 700);
    TGraph *gL = new TGraphErrors(nbins, x_L, y_L, x_error_L, y_error_L);
    gL->Draw("ap");
    TLine *lineL = make_diagonal_line(cL);
    lineL->Draw();
    
    TCanvas *cR = new TCanvas("costheta R", "cR", 80, 80, 700, 700);
    TGraph *gR = new TGraphErrors(nbins, x_R, y_R, x_error_R, y_error_R);
    gR->Draw("ap");
    TLine *lineR = make_diagonal_line(cR);
    lineR->Draw();
*/
}
