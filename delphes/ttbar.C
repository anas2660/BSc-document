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

#include "TLorentzVector.h"

#define GET_BRANCH(x)   ((GenParticle *)branchGenParticle->At(x))
#define GET_PID(x)      (GET_BRANCH(x)->PID)
#define GET_P4(x)       (GET_BRANCH(x)->P4())
#define ABS_PID(x)      (abs(GET_PID(x)))
#define GET_ELECTRON(x) ((Electron *)branchElectron->At(x))
#define GET_MUON(x)     ((Muon *)branchMuon->At(x))

//------------------------------------------------------------------------------

TH1 *histCTStarGen;

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
  //TClonesArray *branchPhoton      = treeReader->UseBranch("Photon");

  // Book histograms
  TH1 *histMissET = new TH1F("missET",  "Missing E_{T}", 60, 0., 300.);
  TH1 *histLeptPT = new TH1F("lept_pt", "lepton P_{T}",  60, 0., 300.);
  TH1 *histJet0PT = new TH1F("jet0_pt", "jet0 P_{T}",    60, 0., 300.);
  TH1 *histJet1PT = new TH1F("jet1_pt", "jet1 P_{T}",    60, 0., 300.);
  TH1 *histJet2PT = new TH1F("jet2_pt", "jet2 P_{T}",    60, 0., 300.);
  TH1 *histJet3PT = new TH1F("jet3_pt", "jet3 P_{T}",    60, 0., 300.);
  TH1 *histNBTag  = new TH1F("NBTag",   "number BTag",   5,  0., 5.);

  histCTStarGen =
      new TH1F("CTStarGen", "generated cos(theta*)", 15, -1., 1.);

  TH1 *histCTStar =
      new TH1F("CTStar", "new cos(theta*)", 15, -1., 1.);

  // Loop over all events
  for (int entry = 0; entry < numberOfEntries; ++entry) {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

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

    if (ctstarGen >= -1.0 && ctstarGen <= 1.0) 
      histCTStarGen->Fill(ctstarGen);
      
    
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

    histLeptPT->Fill(plep.Perp());

    // Missing ET
    MissingET *met = (MissingET *)branchMissingET->At(0);

    histMissET->Fill(met->MET);

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
    histJet0PT->Fill(jet[0]->PT);
    histJet1PT->Fill(jet[1]->PT);
    histJet2PT->Fill(jet[2]->PT);
    histJet3PT->Fill(jet[3]->PT);

    // At least 4 jets with PT > 30 GeV
    if (jet[3]->PT < 30) continue;

    //printf("ETAS: %f, %f, %f, %f... N: %d\n", jet[0]->Eta, jet[1]->Eta, jet[2]->Eta, jet[3]->Eta, branchJet->GetEntries());




    // printf("BTAG: %d, %d, %d, %d\n", jet[0]->BTag, jet[1]->BTag,
    // jet[2]->BTag, jet[3]->BTag );

    //int nBTag = jet[0]->BTag + jet[1]->BTag + jet[2]->BTag + jet[3]->BTag;
    int ngoodbjets = 0;
    int goodbjet[4];
    for (int i = 0; i < 4; i++)
      if (jet[i]->BTag) goodbjet[ngoodbjets++] = i;

    histNBTag->Fill(ngoodbjets);

    if (ngoodbjets != 2) continue;
    //if (met->MET <= 30000.) continue;

    TLorentzVector MeT;
    MeT.SetPtEtaPhiE(met->MET, 0, met->Phi, met->MET);

    // float mtw = sqrt(2.0f * plep.Pt() * MeT.Et() * (1.0f - cos(plep.DeltaPhi(MeT))));
    // if (mtw > 30000.0f) continue;

    // how do we select bjet????????
    TLorentzVector bjet_1 = jet[goodbjet[1]]->P4();

    float m_t = 172.8f;
    float M_W = 80.4f;
    float Meb = (plep + bjet_1).M();
    float costheta  = 2.0f * Meb * Meb / (m_t * m_t - M_W * M_W) - 1.0f;
    // printf("Meb : %f, \tcos theta : %f\n", Meb, costheta);
    if (costheta >= -1.0f && costheta <= 1.0f) {
      histCTStar->Fill(costheta);
      if (hist) hist->Fill(costheta);
    }

  }
    
  // Show histograms
  //if (histsum != NULL) {
  // TODO FIT distributions
  TCanvas *ctruth = new TCanvas("ctruth", "ctruth", 80, 80, 700, 700);


  // TFitResultPtr fit_result_ptr = atlas_costheta->Fit(f, "S");
  
  histCTStarGen->SetLineWidth(3);
  if (!strcmp(inputFile, "./delphes0.root")){
    TF1* f = new TF1("f0", zero_distribution, -1, 1, 1);
    f->SetParNames("c");
    f->SetNpx(1000);
    f->SetLineColor(kRed);
    f->SetLineWidth(2);
    histCTStarGen->Fit(f);
    histCTStarGen->Draw();
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

  if (hist) return;

  TCanvas *c1 = new TCanvas("c1", "c1", 20, 20, 700, 700);
  c1->Divide(2, 2);
  c1->cd(1);
  histJet0PT->Draw();
  c1->cd(2);
  histJet1PT->Draw();
  c1->cd(3);
  histJet2PT->Draw();
  c1->cd(4);
  histJet3PT->Draw();

  TCanvas *c2 = new TCanvas("c2", "c2", 40, 40, 700, 700);
  c2->Divide(2, 2);
  c2->cd(1);
  histMissET->Draw();
  c2->cd(2);
  histLeptPT->Draw();

  TCanvas *c3 = new TCanvas("c3", "c3", 60, 60, 700, 700);
  histNBTag->Draw();

  TCanvas *c4 = new TCanvas("c4", "c4", 80, 80, 700, 700);
  histCTStarGen->Draw();

  TCanvas *c5 = new TCanvas("c5", "c5", 80, 80, 700, 700);
  histCTStar->Draw();
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

  /*
  FCN=0.00831042 FROM MIGRAD    STATUS=CONVERGED     137 CALLS         138 TOTAL
                     EDM=1.61199e-07    STRATEGY= 1  ERROR MATRIX UNCERTAINTY   1.9 per cent
  EXT PARAMETER                                   STEP         FIRST
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
   1  F_0          4.06808e-01   5.35428e+00  -1.56390e-02   3.09175e-04
   2  F_L          5.05306e-01   2.86144e+00   2.19484e-02   2.48891e-04
   3  F_R          8.81551e-02   3.07593e+00  -4.53229e-03   1.91461e-04
                               ERR DEF= 0.5

*/

}

double fitfunc(double* x, double* pars) {
  // Here I just create the histograms I would like to fit with weights.
  // You don't have to do this here, all that is needed is that the histograms
  // are visible in this function.
  /*
  static bool initialized = false;
  static TH1D h1("h1", "", 100, -10, 10);
  static TH1D h2("h2", "", 100, -10, 10);
  if (!initialized) {
    int n = 100000;
    for (int i=0; i<n; ++i) h1.Fill(gRandom->Gaus(-3));
    for (int i=0; i<n; ++i) h2.Fill(gRandom->Gaus(+3));
    h1.Scale(1./n);
    h2.Scale(1./n);
    initialized = true;
  }
  */

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
  //return 4.0*w1*y1/3.0 + 8.0*w2*y2/3.0 + 8.0*w3*y3/3.0;
/*
 FCN=0.00831042 FROM MIGRAD    STATUS=CONVERGED     137 CALLS         138 TOTAL
                     EDM=1.61199e-07    STRATEGY= 1  ERROR MATRIX UNCERTAINTY   1.9 per cent
  EXT PARAMETER                                   STEP         FIRST
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
   1  F_0          4.06808e-01   5.35428e+00  -1.56390e-02   3.09175e-04
   2  F_L          5.05306e-01   2.86144e+00   2.19484e-02   2.48891e-04
   3  F_R          8.81551e-02   3.07593e+00  -4.53229e-03   1.91461e-04
                               ERR DEF= 0.5
*/
}



void costheta(){
#ifdef __CLING__
  gSystem->Load("libDelphes");
#endif


  TH1 *histCTStarGenSum =
      new TH1F("CTStarGenSum", "generated sum cos(theta*)", 15, -1., 1.);

  histCTStarGen0 =
      new TH1F("CTStarGen0", "cos(theta*)", 15, -1., 1.);
  histCTStarGenL =
      new TH1F("CTStarGenL", "cos(theta*)", 15, -1., 1.);
  histCTStarGenR =
      new TH1F("CTStarGenR", "cos(theta*)", 15, -1., 1.);

  ttbar("./delphes0.root", histCTStarGen0);
  ttbar("./delphesL.root", histCTStarGenL);
  ttbar("./delphesR.root", histCTStarGenR);

  histCTStarGen0->Scale(1.0/histCTStarGen0->GetEntries());
  histCTStarGenL->Scale(1.0/histCTStarGenL->GetEntries());
  histCTStarGenR->Scale(1.0/histCTStarGenR->GetEntries());
  histCTStarGen0->SetLineWidth(3);
  histCTStarGenL->SetLineWidth(3);
  histCTStarGenR->SetLineWidth(3);

  TFile* atlas_file     = new TFile("data.root");
  TH1F*  atlas_costheta = atlas_file->Get<TH1F>("hist_costheta");
  int    nbins  = atlas_costheta->GetNbinsX();

  for (int i = 0; i < nbins; i++) {
    histCTStarGen0->SetBinError(i+1, 0);
    histCTStarGenL->SetBinError(i+1, 0);
    histCTStarGenR->SetBinError(i+1, 0);
  }

  atlas_costheta->Scale(1.0/atlas_costheta->GetEntries());

  //atlas_costheta->SetBinContent(atlas_costheta->GetXaxis()->GetNbins()+1, 0);
  //histCTStarGen0->SetBinContent(histCTStarGen0->GetXaxis()->GetNbins()+1, 0);
  //histCTStarGenL->SetBinContent(histCTStarGenL->GetXaxis()->GetNbins()+1, 0);
  //histCTStarGenR->SetBinContent(histCTStarGenR->GetXaxis()->GetNbins()+1, 0);

  TCanvas *c5 = new TCanvas("c5", "c5", 80, 80, 700, 700);
  TF1* f = new TF1("fitfunc", fitfunc, -1, 1, 3);
  f->SetParNames("F_0", "F_L", "F_R");
  f->SetNpx(2500);
  f->SetLineColor(kRed);
  f->SetLineWidth(2);

  //atlas_costheta->Scale(1.0f/128.460461f);

  TFitResultPtr fit_result_ptr = atlas_costheta->Fit(f, "S");
  TFitResult* fit_result = fit_result_ptr.Get();
  printf("CHI2 : %g\n", fit_result->Chi2());

  double parameters[3];
  double f0 = parameters[0] = fit_result->Parameter(0);
  double fl = parameters[1] = fit_result->Parameter(1);
  double fr = parameters[2] = fit_result->Parameter(2);

  const char* filenames[4] = {"dataA.root","dataB.root","dataC.root","dataD.root"};
  TFile* data_files[4];
  TH1F*  data_hist[4];

  double* error_array = (double*) malloc(sizeof(double)*nbins+2);
  double* errors = &error_array[1];
  double* xbar   = (double*)malloc(sizeof(double)*nbins);

  for (int i = 0; i < nbins; i++){
    errors[i] = 0;
    xbar[i]   = atlas_costheta->GetBinContent(i+1);
  }

  for (int i = 0; i < 4; i++) {
    data_files[i] = new TFile(filenames[i]);
    data_hist[i]  = data_files[i]->Get<TH1F>("hist_costheta");
    data_hist[i]->Scale(1.0/data_hist[i]->GetEntries());
    for (int bin = 0; bin < nbins; bin++) {
      double diff = data_hist[i]->GetBinContent(bin+1) - xbar[bin];
      diff = abs(diff);
      diff *= diff;
      errors[bin] += diff;
    }
  }

  for (int i = 0; i < nbins; i++)
    errors[i] = sqrt(errors[i]/4.0);


  /*
  // Draw fit error bars
  double* x     = (double*)malloc(sizeof(double)*nbins); //{1,2,3};
  double* y     = (double*)malloc(sizeof(double)*nbins); //{1,4,9};
  double* err_x = (double*)malloc(sizeof(double)*nbins); //{0,0,0};
  double* err_y = (double*)malloc(sizeof(double)*nbins); //{5,5,5};


  //F_0", "F_L", "F_R
  for (int i = 0; i < nbins; i++) {
    x[i]     =  atlas_costheta->GetBinCenter(i+1);
    y[i]     =  fitfunc(&x[i], parameters);
    err_x[i] =  0;
    err_y[i] =  f0*histCTStarGen0->GetBinError(i);
    err_y[i] += fl*histCTStarGenL->GetBinError(i);
    err_y[i] += fr*histCTStarGenR->GetBinError(i);
  }

   auto errorband = new TGraphErrors(nbins,x,y,err_x,err_y);
  // errorband->SetFillColor(kRed);
  // errorband->SetFillStyle(3010);
  // errorband->Draw("a3 SAME");
  */

  atlas_costheta->SetLineWidth(3);
  atlas_costheta->SetError(error_array);
  atlas_costheta->SetLineColor(kBlue);
  atlas_costheta->Draw("SAME");
  atlas_costheta->ShowBackground()->SetLineColor(kBlue);


  // Draw legend
  auto legend = new TLegend(0.6,0.7125,0.975,0.75);
  legend->AddEntry(atlas_costheta, "ATLAS data cos theta*","l");
  //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  // legend->AddEntry("fitfunc","Delphes fit","l");
  legend->Draw();

  c5->SaveAs("out/delphes_fit.png");

  //TCanvas *cfit = new TCanvas("cfit", "cfit", 80, 80, 700, 700);
//
//  atlas_costheta->Draw();
//
//  errorband->SetFillColor(kRed);
//  errorband->SetFillStyle(3010);
//  errorband->Draw("a3 same");
//
//  legend->Draw();


  TCanvas *c6 = new TCanvas("c6", "c6", 80, 80, 700, 700);
  histCTStarGen0->Draw();
  histCTStarGen0->ShowBackground()->SetLineColor(kBlue);
  c6->SaveAs("out/delphes_ctstar0.png");

  TCanvas *c7 = new TCanvas("c7", "c7", 80, 80, 700, 700);
  histCTStarGenR->Draw(); //"SAME"
  histCTStarGenR->ShowBackground()->SetLineColor(kBlue);
  c7->SaveAs("out/delphes_ctstarR.png");

  TCanvas *c8 = new TCanvas("c8", "c8", 80, 80, 700, 700);
  //histCTStarGenR->SetLineColor(gROOT->GetColor(1)->GetNumber());
  histCTStarGenL->Draw(); //"SAME"
  histCTStarGenL->ShowBackground()->SetLineColor(kBlue);

  c8->SaveAs("out/delphes_ctstarL.png");

  assert(atlas_file);
  assert(atlas_costheta);
  assert(histCTStarGen0);
  assert(histCTStarGenL);
  assert(histCTStarGenR);

}

int main(int argc, char *argv[]) {
  costheta();
  return 0;
}
