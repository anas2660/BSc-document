#include <TFile.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLegend.h>
#include <TGaxis.h>

void efficiency(){
    gStyle->SetOptStat(0);

    TFile* delphes_file = new TFile("delphes.root");
    TH1F* truth_0       = delphes_file->Get<TH1F>("CTStarGen0");
    TH1F* truth_R       = delphes_file->Get<TH1F>("CTStarGenR");
    TH1F* truth_L       = delphes_file->Get<TH1F>("CTStarGenL");
    TH1F* costheta_0    = delphes_file->Get<TH1F>("CTStar0");
    TH1F* costheta_R    = delphes_file->Get<TH1F>("CTStarR");
    TH1F* costheta_L    = delphes_file->Get<TH1F>("CTStarL");
    TH1F* truth_sum     = (TH1F*)truth_0->Clone("Truth_Sum");
    TH1F* costheta_sum  = (TH1F*)costheta_0->Clone("CTStar_Sumold");
    int nbins = truth_0->GetNbinsX();

    TH1F* truth_sum_new               = delphes_file->Get<TH1F>("ctstar_sum_truth");
    TH1F* ctstar_sum                  = delphes_file->Get<TH1F>("ctstar_sum");
    TH1F* ctstar_sum_truth_after_cuts = delphes_file->Get<TH1F>("ctstar_sum_truth_cut");

    truth_0->Sumw2();
    truth_L->Sumw2();
    truth_R->Sumw2();
    truth_sum->Sumw2();
    costheta_0->Sumw2();
    costheta_L->Sumw2();
    costheta_R->Sumw2();    
    costheta_sum->Sumw2();
    ctstar_sum->Sumw2();
    ctstar_sum_truth_after_cuts->Sumw2();
    truth_sum_new->Sumw2();
    
    truth_sum->Add(truth_R);
    truth_sum->Add(truth_L);
    costheta_sum->Add(costheta_R);
    costheta_sum->Add(costheta_L);
    
    costheta_0->Divide(truth_0);
    costheta_L->Divide(truth_L);
    costheta_R->Divide(truth_R);
    //costheta_sum->Divide(truth_sum);
    //ctstar_sum->Divide(ctstar_sum_truth_after_cuts);

    TCanvas *c0 = new TCanvas("costheta0old", "c0", 80, 80, 1400, 500);
    c0->SetTitle("Delphes Scatterplot");
    c0->Divide(2, 2);
    c0->cd(1); costheta_L->Draw("e");
    c0->cd(2); costheta_R->Draw("e");
    c0->cd(3); costheta_0->Draw("e");
    c0->cd(4); costheta_sum->Draw("e");

    c0->SaveAs("../figures/Efficiencyold.png");


    TCanvas *c1 = new TCanvas("costheta0_efficiency", "c1", 80, 80, 700, 500);
    c1->SetTitle("Delphes Efficiency");
    auto a = (TH1F*)ctstar_sum_truth_after_cuts->Clone("a");
    a->Divide(truth_sum_new);
    a->SetTitle("Efficiency of Cuts Per Bin");
    a->GetXaxis()->SetTitle("cos #theta*");
    a->GetYaxis()->SetTitle("Efficiency of Cuts");

    // axis titles
    a->Draw("e");
    a->ShowBackground()->SetLineColor(kBlue);
    // ctstar_sum->Draw("e");

    c1->SaveAs("../figures/Efficiency.png");

    TCanvas *c2 = new TCanvas("costheta0", "c2", 80, 80, 700, 500);
    c2->Divide(2);
    // TITLES
    // Axis titles
    c2->SetTitle("Delphes ");
    c2->cd(1); ctstar_sum->Draw("e");
    c2->cd(2); ctstar_sum_truth_after_cuts->Draw("e");

    //c1->SaveAs("../figures/Efficiency.png");


    //TLegend* leg = new TLegend(0.1,0.7,0.3,0.9);
    ////leg->SetHeader("cos theta* 0");
    //leg->AddEntry(gL,    "costheta* L", "lep");
    //leg->AddEntry(g0,    "costheta* 0", "lep");
    //leg->AddEntry(gR,    "costheta* R", "lep");
    //leg->AddEntry(line0, "Expectation","l r");
    //leg->Draw();


}
