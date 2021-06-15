/*
Simple macro showing how to access branches from the delphes output root file,
loop over events, and plot simple quantities such as the jet pt and the
di-electron invariant mass.

root -l examples/Example1.C'("delphes_output.root")'
*/

#include <cassert>
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#endif
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

#include <TH1F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TClonesArray.h>
#include <TH2F.h>


#include "TLorentzVector.h"

#define GET_BRANCH(x)   ((GenParticle *)branchGenParticle->At(x))
#define GET_PID(x)      (GET_BRANCH(x)->PID)
#define GET_P4(x)       (GET_BRANCH(x)->P4())
#define ABS_PID(x)      (abs(GET_PID(x)))
#define GET_ELECTRON(x) ((Electron *)branchElectron->At(x))
#define GET_MUON(x)     ((Muon *)branchMuon->At(x))

//------------------------------------------------------------------------------

TH1   *histCTStarGen;
TFile *out_file;

TH2D  *ctRecVsTrue;

TH1F  *ctstar_sum;
TH1F  *ctstar_sum_truth;
TH1F  *ctstar_sum_truth_after_cuts;


///////////////////////////////////////////////////////////////////////////////
//                               Distributions                               //
///////////////////////////////////////////////////////////////////////////////

double left_distribution(double *x, double *par) {
  double y;
  y = 1 - x[0];
  y = y * y * 3.0 / 8.0;
  return y*par[0];
}

double zero_distribution(double *x, double *par) {
  double y;
  y = x[0];
  y = (1 - y*y) * 3.0 / 4.0;
  return y*par[0];
}

double right_distribution(double *x, double *par) {
  double y;
  y = 1 + x[0];
  y = y * y * 3.0 / 8.0;
  return y*par[0];
}

void ttbar(const char *inputFile, TH1 *hist) {
  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries     = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchEvent       = treeReader->UseBranch("Event");
  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchMissingET   = treeReader->UseBranch("MissingET");
  TClonesArray *branchJet         = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron    = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon        = treeReader->UseBranch("Muon");
  TClonesArray *branchWeight      = treeReader->UseBranch("HepMCEvent");
  //TClonesArray *branchPhoton      = treeReader->UseBranch("Photon");

  // Book histograms
  // TH1 *histMissET = new TH1F("missET",  "Missing E_{T}", 60, 0., 300.);
  // TH1 *histLeptPT = new TH1F("lept_pt", "lepton P_{T}",  60, 0., 300.);
  // TH1 *histJet0PT = new TH1F("jet0_pt", "jet0 P_{T}",    60, 0., 300.);
  // TH1 *histJet1PT = new TH1F("jet1_pt", "jet1 P_{T}",    60, 0., 300.);
  // TH1 *histJet2PT = new TH1F("jet2_pt", "jet2 P_{T}",    60, 0., 300.);
  // TH1 *histJet3PT = new TH1F("jet3_pt", "jet3 P_{T}",    60, 0., 300.);
  // TH1 *histNBTag  = new TH1F("NBTag",   "number BTag",   5,  0., 5.);

  // histCTStarGen=
  //     new TH1F("CTStarGen", "truth cos(theta*)", 15, -1., 1.);


  if (!strcmp(inputFile, "./delphes0.root"))
    histCTStarGen =
      new TH1F("CTStarGen0", "Predicted #Gamma_{0} Fitted to Generated cos(theta*)", 15, -1., 1.);
  if (!strcmp(inputFile, "./delphesL.root"))
    histCTStarGen =
      new TH1F("CTStarGenL", "Predicted #Gamma_{L} Fitted to Generated cos(theta*)", 15, -1., 1.);
  if (!strcmp(inputFile, "./delphesR.root"))
    histCTStarGen =
      new TH1F("CTStarGenR", "Predicted #Gamma_{R} Fitted to Generated cos(theta*)", 15, -1., 1.);
  histCTStarGen->GetXaxis()->SetTitle("cos #theta*");
  histCTStarGen->GetYaxis()->SetTitle("Events");
  //histCTStarGen->GetYaxis()->SetTitleOffset(1.);

  TH1 *histCTStar =
      new TH1F("CTStar", "new cos(theta*)", 15, -1., 1.);

  // Loop over all events
  for (int entry = 0; entry < numberOfEntries; ++entry) {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);


    // Load weight for event
    //// {
    ////   int count = branchWeight->GetEntries();
    ////   if (count) printf("Weight count: %d\n", count);
    //// }
    //printf("Weight count: %d\n", branchWeight->GetEntries());
    //float weight = ((Weight *)branchWeight->At(0))->Weight;
    //float weight = ((Weight *)branchWeight->At(entry))->Weight;
    //printf("Weight: %g\n", weight);

    // Treat truth information:
    // Cumbersome, since daughter pointers seem no to be set correctly
    int ibm  = -1;
    int ibp  = -1;
    int ilep = -1;
    int npar = branchGenParticle->GetEntries();

    for (int ipar = 0; ipar < npar; ++ipar) {
      // Find lepton (e or mu) from W decay
      auto ipar_branch = GET_BRANCH(ipar);
      auto PID = ABS_PID(ipar);
      if (PID == 11 || PID == 13) {
        int M1 = ipar_branch->M1;
        if (ABS_PID(M1) == 24) ilep = ipar;
      }

      // Find b from t decay
      if (PID == 5) { // b particle
        int M1 = (ipar_branch)->M1;
        if (M1 < 0) continue;
        if (ABS_PID(M1) == 6) { // top particle
          if (GET_PID(ipar) < 0) ibm = ipar;
          else ibp = ipar;
        }
      }

      // When found what we came for, break out to save time
      if (ibm >= 0 && ibp >= 0 && ilep >= 0) break;
    }

    // No lepton from W found!!! This should not happen. Skip event if it does!
    if (ilep < 0) break;

    // Select b or anti-b depending on charge of lepton
    int ib = ibm;
    if (GET_PID(ilep) < 0) ib = ibp;

    // Calculate cos(theta*)
    double Meb2      = (GET_P4(ilep) + GET_P4(ib)).M2();
    double ctstarGen = 2. * Meb2 / (172.5 * 172.5 - 80.1 * 80.1) - 1.;

    if (ctstarGen >= -1.0 && ctstarGen <= 1.0){
      histCTStarGen->Fill(ctstarGen);
      ctstar_sum_truth->Fill(ctstarGen);
    }

    
    //if (hist) hist->Fill(ctstarGen);

    // Turn to simulated events...
    int nelec   = branchElectron->GetEntries();
    int nGoodEl = 0;
    int good_electron_id = -1;
    for (int ielec = 0; ielec < nelec; ++ielec) {
      Electron *elec = GET_ELECTRON(ielec);
      float     eta  = elec->Eta;
      if (elec->PT > 30. && fabsf(eta) < 2.47 && (fabsf(eta) < 1.37 || fabsf(eta) > 1.52)) {
        nGoodEl++;
        if(good_electron_id == -1) good_electron_id = ielec;
      }
    }

    int nmuon        = branchMuon->GetEntries();
    int nGoodMu      = 0;
    int good_muon_id = -1;

    for (int imuon = 0; imuon < nmuon; ++imuon) {
      Muon *muon = GET_MUON(imuon);
      if (muon->PT > 30. && fabsf(muon->Eta) < 2.5){
        nGoodMu++;
        if(good_muon_id == -1) good_muon_id = imuon;
      }
    }


    // Exactly 1 lepton with pT > 30 GeV
    if (nGoodEl + nGoodMu != 1) continue;

    TLorentzVector plep;
    if (nGoodEl == 1) plep = GET_ELECTRON(good_electron_id)->P4();
    if (nGoodMu == 1) plep = GET_MUON(good_muon_id)->P4();

    ///////histLeptPT->Fill(plep.Perp());

    // Missing ET
    MissingET *met = (MissingET *)branchMissingET->At(0);

    //////histMissET->Fill(met->MET);

    // AT least 30 GeV misisng ET
    if (met->MET < 30.) continue;

    // At least 4 jets
    if (branchJet->GetEntries() < 4) continue;

    // Take first four jet (they are ordered by pT)
    Jet *jet[4];
    int jet_count = 0;
    for (int i = 0; i < branchJet->GetEntries(); i++) {
      if (jet_count == 4) break;
      Jet* current_jet =  (Jet*)branchJet->At(i);
      if (fabsf(current_jet->Eta) < 2.5)
        jet[jet_count++] = current_jet;
    }

    if (jet_count != 4) continue;

    // printf("a1: 0x%x, sizeof: %lu, a2: 0x%x\n", branchJet->At(0), sizeof(Jet), branchJet->At(1));


    // Plot jet transverse momentum
    /////// histJet0PT->Fill(jet[0]->PT);
    /////// histJet1PT->Fill(jet[1]->PT);
    /////// histJet2PT->Fill(jet[2]->PT);
    /////// histJet3PT->Fill(jet[3]->PT);

    // At least 4 jets with PT > 30 GeV
    if (jet[3]->PT < 30) continue;

    //int nBTag = jet[0]->BTag + jet[1]->BTag + jet[2]->BTag + jet[3]->BTag;
    int ngoodbjets = 0;
    int goodbjet[4];
    for (int i = 0; i < 4; i++)
      if (jet[i]->BTag) goodbjet[ngoodbjets++] = i;

    ///// histNBTag->Fill(ngoodbjets);

    if (ngoodbjets != 2) continue;
    //if (met->MET <= 30000.) continue;

    ///////TLorentzVector MeT;
    ///////MeT.SetPtEtaPhiE(met->MET, 0, met->Phi, met->MET);

    // float mtw = sqrt(2.0f * plep.Pt() * MeT.Et() * (1.0f - cos(plep.DeltaPhi(MeT))));
    // if (mtw > 30000.0f) continue;

    TLorentzVector bjet_1 = jet[goodbjet[1]]->P4();

    // puts("========"); // (GenParticle *)
    // printf("a: %d, b: %d, c: %d\n", jet[goodbjet[1]]->BTag, jet[goodbjet[1]]->BTagAlgo, jet[goodbjet[1]]->BTagPhys);

    // printf("COUNT: %d\n", jet[goodbjet[1]]->Particles.GetEntries());
    // for (int i = 0; i < jet[goodbjet[1]]->Particles.GetEntries(); i++) {
    //   int type =  ((GenParticle*)jet[goodbjet[1]]->Particles.At(i))->PID;
    //   printf("type: %d\n", type);
    // }

    // printf("COUNT2: %d\n", jet[goodbjet[0]]->Particles.GetEntries());
    // for (int i = 0; i < jet[goodbjet[0]]->Particles.GetEntries(); i++) {
    //   int type =  ((GenParticle*)jet[goodbjet[0]]->Particles.At(i))->PID;
    //   printf("type2: %d\n", type);
    // }

    // printf("COUNT3: %d\n", jet[goodbjet[1]]->Constituents.GetEntries());
    // for (int i = 0; i < jet[goodbjet[1]]->Constituents.GetEntries(); i++) {
    //   int type =  ((GenParticle*)jet[goodbjet[1]]->Constituents.At(i))->PID;
    //   printf("type3: %d\n", type);
    // }

    // puts(jet[goodbjet[0]]->Constituents.At(0)->GetName());
    // printf("COUNT4: %d\n", jet[goodbjet[0]]->Constituents.GetEntries());
    // for (int i = 0; i < jet[goodbjet[0]]->Constituents.GetEntries(); i++) {
    //   int type =  ((GenParticle*)jet[goodbjet[0]]->Constituents.At(i))->PID;
    //   printf("type4: %d\n", type);
    // }

    //puts("========");

    double m_t = 172.8f;
    double M_W = 80.4f;
    double Meb = (plep + bjet_1).M();
    double costheta  = 2.0f * Meb * Meb / (m_t * m_t - M_W * M_W) - 1.0f;
    // printf("Meb : %f, \tcos theta : %f\n", Meb, costheta);
    if (costheta >= -1.0f && costheta <= 1.0f) {
      histCTStar->Fill(costheta);
      ctRecVsTrue->Fill(ctstarGen, costheta);
      ctstar_sum->Fill(costheta);
      ctstar_sum_truth_after_cuts->Fill(ctstarGen);
      if (hist) hist->Fill(costheta);
    }
  }
    
  // Show histograms
  //if (histsum != NULL) {
  // TODO FIT distributions
  TCanvas *ctruth = new TCanvas("ctruth", "ctruth", 80, 80, 900, 700);


  // TFitResultPtr fit_result_ptr = atlas_costheta->Fit(f, "S");
  
  histCTStarGen->SetLineWidth(3);
  if (!strcmp(inputFile, "./delphes0.root")){
    TF1* f = new TF1("f0", zero_distribution, -1, 1, 1);
    f->SetParNames("c");
    f->SetNpx(1000);
    f->SetLineColor(kRed);
    f->SetLineWidth(2);
    histCTStarGen->Fit(f);
    histCTStarGen->GetYaxis()->SetRangeUser(0, 2100);
    histCTStarGen->Draw();
    histCTStarGen->Write();
    ctruth->SaveAs("out/delphes_gen0.png");
  }

  //histCTStarGen->Scale(0.687);
  if (!strcmp(inputFile, "./delphesL.root")){
    TF1* f = new TF1("fl", left_distribution, -1, 1, 1);
    f->SetParNames("c");
    f->SetNpx(1000);
    f->SetLineColor(kRed);
    f->SetLineWidth(2);
    histCTStarGen->Fit(f);
    histCTStarGen->Draw();
    histCTStarGen->Write();
    ctruth->SaveAs("out/delphes_genL.png");
  }
  //histCTStarGen->Scale(0.311);
  if (!strcmp(inputFile, "./delphesR.root")){
    TF1* f = new TF1("fr", right_distribution, -1, 1, 1);
    f->SetParNames("c");
    f->SetNpx(1000);
    f->SetLineColor(kRed);
    f->SetLineWidth(2);
    histCTStarGen->Fit(f);
    histCTStarGen->Draw();
    histCTStarGen->Write();
    ctruth->SaveAs("out/delphes_genR.png");
  }

  //histCTStarGen->Scale(0.0017);
  //  histsum->Add(histCTStarGen);
  //  return;
  //}
  ctruth->Close();
  delete ctruth;

  //TCanvas *c4 = new TCanvas("c4", "c4", 80, 80, 700, 700);
  //histCTStarGen->Draw();

  //if (hist) return;

  // TCanvas *c1 = new TCanvas("c1", "c1", 20, 20, 700, 700);
  // c1->Divide(2, 2);
  // c1->cd(1);
  // histJet0PT->Draw();
  // c1->cd(2);
  // histJet1PT->Draw();
  // c1->cd(3);
  // histJet2PT->Draw();
  // c1->cd(4);
  // histJet3PT->Draw();

  // TCanvas *c2 = new TCanvas("c2", "c2", 40, 40, 700, 700);
  // c2->Divide(2, 2);
  // c2->cd(1);
  // histMissET->Draw();
  // c2->cd(2);
  // histLeptPT->Draw();

  // TCanvas *c3 = new TCanvas("c3", "c3", 60, 60, 700, 700);
  // histNBTag->Draw();

  // TCanvas *c4 = new TCanvas("c4", "c4", 80, 80, 700, 700);
  // histCTStarGen->Draw();

  // TCanvas *c5 = new TCanvas("c5", "c5", 80, 80, 700, 700);
  // histCTStar->Draw();
}


TH1 *histCTStarGen0;
TH1 *histCTStarGenL;
TH1 *histCTStarGenR;

// define a function with 3 parameters
double fitf(double *x, double *par) {
  double delphes_R, delphes_0, delphes_L;
  delphes_0 = histCTStarGen0->Interpolate(x[0]);
  delphes_L = histCTStarGenL->Interpolate(x[0]);
  delphes_R = histCTStarGenR->Interpolate(x[0]);
  return delphes_0*par[0] + delphes_L*par[1] + delphes_R*par[2];

}

double fitfunc(double* x, double* pars) {
  const double xx = x[0]; // use x[1] to get 2nd dimension, x[2] for 3rd ...
  // the fit parameters, i.e. the histogram weights
  const double w1 = pars[0];
  const double w2 = pars[1];
  const double w3 = pars[2];

  // get content of the histograms for this point
  const double y1 = histCTStarGen0->GetBinContent(histCTStarGen0->GetXaxis()->FindFixBin(xx));
  const double y2 = histCTStarGenL->GetBinContent(histCTStarGenL->GetXaxis()->FindFixBin(xx));
  const double y3 = histCTStarGenR->GetBinContent(histCTStarGenR->GetXaxis()->FindFixBin(xx));

  return w1*y1 + w2*y2 + w3*y3;
}



void costheta(){
#ifdef __CLING__
  gSystem->Load("libDelphes");
#endif

  gStyle->SetOptStat(0);
  ctRecVsTrue                 = new TH2D("cthRecVsTrue", "Reconstructed cos(#theta*) vs Generated cos(#theta*)", 15, -1., 1., 15, -1., 1.);
  ctstar_sum                  = new TH1F("ctstar_sum", "sum cos(theta*)", 15, -1, 1);
  ctstar_sum_truth            = new TH1F("ctstar_sum_truth", "sum cos(theta*) truth", 15, -1, 1);
  ctstar_sum_truth_after_cuts = new TH1F("ctstar_sum_truth_cut", "sum truth cos(theta*) cut", 15, -1, 1);;

  ctRecVsTrue->GetXaxis()->SetTitle("Generated cos #theta*");
  ctRecVsTrue->GetYaxis()->SetTitle("Reconstructed cos #theta*");

  TFile* atlas_file     = new TFile("data.root");
  TH1F*  atlas_costheta = atlas_file->Get<TH1F>("hist_costheta");
  int    nbins  = atlas_costheta->GetNbinsX();

  TH1 *histCTStarGenSum =
      new TH1F("CTStarGenSum", "generated sum cos(theta*)", 15, -1., 1.);

  histCTStarGen0 = new TH1F("CTStar0", "cos(theta*)", 15, -1., 1.);
  histCTStarGenL = new TH1F("CTStarL", "cos(theta*)", 15, -1., 1.);
  histCTStarGenR = new TH1F("CTStarR", "cos(theta*)", 15, -1., 1.);

  out_file = new TFile("delphes.root", "recreate");
  ttbar("./delphes0.root", histCTStarGen0);
  ttbar("./delphesL.root", histCTStarGenL);
  ttbar("./delphesR.root", histCTStarGenR);
  histCTStarGen0->Write();
  histCTStarGenL->Write();
  histCTStarGenR->Write();
  ctstar_sum->Write();
  ctstar_sum_truth->Write();
  ctstar_sum_truth_after_cuts->Write();
      

  histCTStarGen0->Scale(1.0/histCTStarGen0->GetEntries());
  histCTStarGenL->Scale(1.0/histCTStarGenL->GetEntries());
  histCTStarGenR->Scale(1.0/histCTStarGenR->GetEntries());
  histCTStarGen0->SetLineWidth(3);
  histCTStarGenL->SetLineWidth(3);
  histCTStarGenR->SetLineWidth(3);

  double* error_array = (double*) malloc(sizeof(double)*nbins+2);
  double* errors = &error_array[1];
  for (int i = 0; i < nbins; i++)
    errors[i] = sqrt(atlas_costheta->GetBinContent(i+1))/atlas_costheta->GetEntries();
  atlas_costheta->Scale(1.0/atlas_costheta->GetEntries());

  TCanvas *c5 = new TCanvas("c5", "c5", 80, 80, 700, 700);
  TF1* f = new TF1("fitfunc", fitfunc, -1, 1, 3);
  f->SetParNames("F_0", "F_L", "F_R");
  f->SetNpx(2500);
  f->SetLineColor(kRed);
  f->SetLineWidth(2);

  //atlas_costheta->Scale(1.0f/128.460461f);

  TFitResultPtr fit_result_ptr = atlas_costheta->Fit(f, "S");
  TFitResult*   fit_result     = fit_result_ptr.Get();
  printf("CHI2 : %g\n", fit_result->Chi2());

  double parameters[3];
  double f0 = parameters[0] = fit_result->Parameter(0);
  double fl = parameters[1] = fit_result->Parameter(1);
  double fr = parameters[2] = fit_result->Parameter(2);

  atlas_costheta->SetError(error_array);
  atlas_costheta->SetLineWidth(3);
  atlas_costheta->SetLineColor(kBlue);
  atlas_costheta->Draw("e SAME");
  atlas_costheta->ShowBackground()->SetLineColor(kBlue);

  atlas_costheta->Write();
  out_file->Close();
  delete out_file;

  // Draw legend
  auto legend = new TLegend(0.6,0.7125,0.975,0.75);
  legend->AddEntry(atlas_costheta, "ATLAS data cos theta*","l");
  //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  // legend->AddEntry("fitfunc","Delphes fit","l");
  legend->Draw();


  double* hist_fit_data  = (double*)malloc(sizeof(double)*nbins+2);
  double* hist_fit_error = (double*)malloc(sizeof(double)*nbins+2);

  for (int i = 0; i < nbins+2; i++) {
    hist_fit_data[i]  = f0 * histCTStarGen0->GetBinContent(i);
    hist_fit_data[i] += fl * histCTStarGenL->GetBinContent(i);
    hist_fit_data[i] += fr * histCTStarGenR->GetBinContent(i);

    // Through error propagation sqrt(e3^2 * x3^2 + e2^2 * x2^2 + e1^2 * x1^2)
    double e1 = histCTStarGen0->GetBinError(i);
    double e2 = histCTStarGenL->GetBinError(i);
    double e3 = histCTStarGenR->GetBinError(i);
    double error = sqrt(e3*e3 * fr*fr + e2*e2 * fl*fl + e1*e1 * f0*f0);
    hist_fit_error[i] = error;
  }
  
  TCanvas* cfit  = new TCanvas("cfit", "cfit", 80, 80, 1200, 700);
  TH1F* hist_fit = new TH1F("hist_name", "Delphes cos#theta* Fit to Atlas Data", nbins, -1, 1);
  hist_fit->SetContent(hist_fit_data);
  hist_fit->SetError(hist_fit_error);

  hist_fit->GetXaxis()->SetTitle("cos#theta*");
  hist_fit->GetYaxis()->SetRangeUser(0, 0.13);

  hist_fit->SetLineWidth(3);
  hist_fit->SetLineColor(kRed);
  hist_fit->Draw("E1");

  atlas_costheta->SetTitle("cos#theta*");
  atlas_costheta->Draw("E1 SAME");
  atlas_costheta->ShowBackground()->SetLineColor(kBlue);
  //atlas_costheta->get

    //TLegend* legend2 = new TLegend(0.75,0.65,0.975,0.75);
  TLegend* legend2 = new TLegend(0.75, 0.775, 0.98, 0.935);
  legend2->AddEntry(atlas_costheta, "ATLAS data cos#theta*",     "l");
  legend2->AddEntry(hist_fit,       "Fitted Delphes cos#theta*", "l");
  // legend2->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  // legend2->AddEntry("fitfunc","Delphes fit","l");
  legend2->Draw();

  cfit->SaveAs("out/delphes_fit.png");

  
  TCanvas *c6 = new TCanvas("c6", "c6", 80, 80, 900, 650);
  histCTStarGen0->GetXaxis()->SetTitle("cos #theta*");
  histCTStarGen0->GetYaxis()->SetTitle("Events");
  histCTStarGen0->SetTitle("");
  //histCTStarGen0->SetTitleSize(0);
  histCTStarGen0->GetYaxis()->SetRangeUser(0, 0.125);
  histCTStarGen0->Draw("E1");
  histCTStarGen0->ShowBackground()->SetLineColor(kBlue);
  c6->SaveAs("out/delphes_ctstar0.png");

  TCanvas *c7 = new TCanvas("c7", "c7", 80, 80, 900, 650);
  histCTStarGenR->GetYaxis()->SetTitle("Events");
  histCTStarGenR->SetTitle("");
  //histCTStarGenR->SetTitleSize(0);
  histCTStarGenR->GetYaxis()->SetRangeUser(0, 0.12);
  histCTStarGenR->GetXaxis()->SetTitle("cos #theta*");
  histCTStarGenR->Draw("E1"); //"SAME"
  histCTStarGenR->ShowBackground()->SetLineColor(kBlue);
  c7->SaveAs("out/delphes_ctstarR.png");

  TCanvas *c8 = new TCanvas("c8", "c8", 80, 80, 900, 650);
  histCTStarGenL->GetXaxis()->SetTitle("cos #theta*");
  histCTStarGenL->GetYaxis()->SetTitle("Events");
  histCTStarGenL->SetTitle("");
  //histCTStarGenL->SetTitleSize(0);
  histCTStarGenL->GetYaxis()->SetRangeUser(0, 0.17);
  histCTStarGenL->Draw("E1"); //"SAME"
  histCTStarGenL->ShowBackground()->SetLineColor(kBlue);
  c8->SaveAs("out/delphes_ctstarL.png");

  gStyle->SetPalette(kRainBow);
  TCanvas * c2d = new TCanvas("c2d","c2d",40,40,900,800);
  ctRecVsTrue->GetYaxis()->SetTitleOffset(1.2f);
  ctRecVsTrue->Draw("colz");

  //ctRecVsTrue->Draw("cont4z");
}

int main(int argc, char *argv[]) {
  costheta();
  return 0;
}
