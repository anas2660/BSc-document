/*
Simple macro showing how to access branches from the delphes output root file,
loop over events, and plot simple quantities such as the jet pt and the
di-electron invariant mass.

root -l examples/Example1.C'("delphes_output.root")'
*/

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
  TH1 *histMissET = new TH1F("missET", "Missing E_{T}", 60, 0., 300.);
  TH1 *histLeptPT = new TH1F("lept_pt", "lepton P_{T}", 60, 0., 300.);
  TH1 *histJet0PT = new TH1F("jet0_pt", "jet0 P_{T}", 60, 0., 300.);
  TH1 *histJet1PT = new TH1F("jet1_pt", "jet1 P_{T}", 60, 0., 300.);
  TH1 *histJet2PT = new TH1F("jet2_pt", "jet2 P_{T}", 60, 0., 300.);
  TH1 *histJet3PT = new TH1F("jet3_pt", "jet3 P_{T}", 60, 0., 300.);
  TH1 *histNBTag  = new TH1F("NBTag", "number BTag", 5, 0., 5.);

  TH1 *histCTStarGen =
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
      if (PID == 5 ) { // b particle
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

    histCTStarGen->Fill(ctstarGen);
    //if (hist) hist->Fill(ctstarGen);

    // Turn to simulated events...
    int nelec   = branchElectron->GetEntries();
    int nGoodEl = 0;
    for (int ielec = 0; ielec < nelec; ++ielec) {
      Electron *elec = GET_ELECTRON(ielec);
      if (elec->PT > 30.) nGoodEl++;
    }

    int nmuon = branchMuon->GetEntries();
    int nGoodMu = 0;
    for (int imuon = 0; imuon < nmuon; ++imuon) {
      Muon *muon = GET_MUON(imuon);
      if (muon->PT > 30.) nGoodMu++;
    }

    // Exactly 1 lepton with pT > 30 GeV
    if (nGoodEl + nGoodMu != 1) continue;

    TLorentzVector plep;
    if (nGoodEl == 1) plep = GET_ELECTRON(0)->P4();
    if (nGoodMu == 1) plep = GET_MUON(0)->P4();

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
    for (int i = 0; i < 4; i++)
      jet[i] = (Jet *)branchJet->At(i);

    // printf("a1: 0x%x, sizeof: %lu, a2: 0x%x\n", branchJet->At(0), sizeof(Jet), branchJet->At(1));


    // Plot jet transverse momentum
    histJet0PT->Fill(jet[0]->PT);
    histJet1PT->Fill(jet[1]->PT);
    histJet2PT->Fill(jet[2]->PT);
    histJet3PT->Fill(jet[3]->PT);

    // At least 4 jets with PT > 30 GeV
    if (jet[3]->PT < 30) continue;

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
    histCTStar->Fill(costheta);
    if (hist) hist->Fill(costheta);

  }
    
  // Show histograms
  //if (histsum != NULL) {
  //  if (!strcmp(inputFile, "./delphes0.root"))
  //    histCTStarGen->Scale(0.687);
  //  if (!strcmp(inputFile, "./delphesL.root"))
  //    histCTStarGen->Scale(0.311);
  //  if (!strcmp(inputFile, "./delphesR.root"))
  //    histCTStarGen->Scale(0.0017);
  //  histsum->Add(histCTStarGen);
  //  return;
  //}

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

  return 4.0*w1*y1/3.0 + 8.0*w2*y2/3.0 + 8.0*w3*y3/3.0;
}



void costheta(){
#ifdef __CLING__
  gSystem->Load("libDelphes");
#endif


  TH1 *histCTStarGenSum =
      new TH1F("CTStarGenSum", "generated sum cos(theta*)", 15, -1., 1.);

  histCTStarGen0 =
      new TH1F("CTStarGen0", "generated cos(theta*)", 15, -1., 1.);
  histCTStarGenL =
      new TH1F("CTStarGenL", "generated cos(theta*)", 15, -1., 1.);
  histCTStarGenR =
      new TH1F("CTStarGenR", "generated cos(theta*)", 15, -1., 1.);

  ttbar("./delphes0.root", histCTStarGen0);
  ttbar("./delphesL.root", histCTStarGenL);
  ttbar("./delphesR.root", histCTStarGenR);

  TFile* atlas_file = new TFile("data.root");
  TH1F* atlas_costheta = atlas_file->Get<TH1F>("hist_costheta");
  atlas_costheta->SetBinContent(atlas_costheta->GetXaxis()->GetNbins()+1, 0);

  histCTStarGen0->SetBinContent(histCTStarGen0->GetXaxis()->GetNbins()+1, 0);
  histCTStarGenL->SetBinContent(histCTStarGenL->GetXaxis()->GetNbins()+1, 0);
  histCTStarGenR->SetBinContent(histCTStarGenR->GetXaxis()->GetNbins()+1, 0);

  TCanvas *c5 = new TCanvas("c5", "c5", 80, 80, 700, 700);
  TF1* f = new TF1("fitfunc", fitfunc, -1, 1, 3);
  f->SetParNames("F_0", "F_L", "F_R");

  //atlas_costheta->Scale(1.0f/128.460461f);

  atlas_costheta->Fit(f, "L");
  atlas_costheta->Draw();
  TCanvas *c6 = new TCanvas("c6", "c6", 80, 80, 700, 700);
  histCTStarGen0->Draw();
  //c6->Divide(3);
  //c6->cd(0);
  ///c6->cd(1);
  ///histCTStarGenL->Draw();
  ///c6->cd(2);
  ///histCTStarGenR->Draw();
  //f->SetParameters(0.687, 0.311, 0.0017);
//f->Draw();

  //histCTStarGen0->Scale(0.687);
  //histCTStarGenL->Scale(0.311);
  //histCTStarGenR->Scale(0.0017);
//
//  histCTStarGenSum->Add(histCTStarGen0);
//  histCTStarGenSum->Add(histCTStarGenL);
//  histCTStarGenSum->Add(histCTStarGenR);
//  histCTStarGenSum->Draw();


  //TF1 *func = new TF1("fit", fitf, -1, 1);
  //func->SetParameters(1.0, 1.0, 1.0);
  //func->SetParNames("F_0", "F_L", "F_R");
//
//  TFitResultPtr resultptr = atlas_costheta->Fit("fit");
//  //TFitResult* result = resultptr.Get();

  //TF1* f = new TF1("fitfunc", fitfunc, -1, 1, 3);
  //atlas_costheta->Fit(f, "L");

  //TCanvas *c5 = new TCanvas("c5", "c5", 80, 80, 700, 700);
  //histCTStarGenSum->Add(histCTStarGen0);
  //histCTStarGenSum->Add(histCTStarGenL);
  //histCTStarGenSum->Add(histCTStarGenR);
  //histCTStarGenSum->Draw();
}

int main(int argc, char *argv[]) {
  costheta();
  return 0;
}
