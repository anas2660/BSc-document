/////////////////////////////////////////////////////////////
//// TTbarAnalysis code
//// Author: ATLAS Collaboration (2019)
////
////
//// DISCLAIMER:
//// This Software is intended for educational use only!
//// Under no circumstances does it qualify to reproduce actual ATLAS analysis
/// results or produce publishable results!
/////////////////////////////////////////////////////////////

#include <TMathBase.h>
#include <cstdio>
#include <math.h>
#define TTbarAnalysis_cxx

#include "TTbarAnalysis.h"
#include "TTbarAnalysisHistograms.h"
#include <cstring>
#include <iostream>
#include <string>

#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TStyle.h>

#define electron_type 11
#define muon_type     13

string name;

void TTbarAnalysis::Begin(TTree *) {
  nEvents = 0; 
  N1 = N3 = 0;
  NTotal = 0;
}

void TTbarAnalysis::SlaveBegin(TTree *) {
  TString option = GetOption();
  printf("Starting analysis with process option: %s \n", option.Data());

  name = option;

  define_histograms();

  FillOutputList();
}

Bool_t TTbarAnalysis::Process(Long64_t entry) {
  fChain->GetTree()->GetEntry(entry);
  nEvents++;
  if (nEvents % 50000 == 0)
    std::cout << "Analysed a total of: " << nEvents << " events out of " << fChain->GetTree()->GetEntries() << " in this sample" << std::endl;

  if (fChain->GetTree()->GetEntries() > 0) {
    // **********************************************************************************************//
    // Begin analysis selection, largely based on: ATLAS Collaboration, JHEP 11
    // (2017) 191           //
    // **********************************************************************************************//

    // Scale factors (adding b-tagging as it is used)
    Float_t scaleFactor = scaleFactor_ELE * scaleFactor_MUON * scaleFactor_LepTRIGGER * scaleFactor_PILEUP * scaleFactor_BTAG;

    // MC weight
    Float_t m_mcWeight = mcWeight;

    // read input option
    TString option = GetOption();
    if (option.Contains("single")) {
      m_mcWeight = (mcWeight / TMath::Abs(mcWeight));
    } // set to 1 or -1 for single top MCs

    // Total weight
    Float_t weight   = scaleFactor * m_mcWeight;
    float   costheta = -2;

    // Make difference between data and MC
    if (option.Contains("data")) {
      weight = 1.;
    }

    // cut on at least 4 jets
    if (jet_n > 3) {

      // MET > 30 GeV
      if (met_et > 30000.) {

        // Preselection cut for electron/muon trigger, Good Run List, and good
        // vertex
        if (trigE || trigM) {

          // Preselection of good leptons
          int goodlep_index = 0;
          int goodlep_n     = 0;
          int lep_index     = 0;

          for (unsigned int i = 0; i < lep_n; i++) {
            TLorentzVector leptemp;
            leptemp.SetPtEtaPhiE(lep_pt->at(i) / 1000., lep_eta->at(i), lep_phi->at(i), lep_E->at(i) / 1000.);

            // Lepton is Tight
            if (lep_isTightID->at(i)) {
              // Lepton is isolated and hard pT
              if (lep_pt->at(i) > 30000. && ((lep_ptcone30->at(i) / lep_pt->at(i)) < 0.15) && ((lep_etcone20->at(i) / lep_pt->at(i)) < 0.15)) {
                // electron selection in fiducial region excluding candidates in
                // the transition region between the barrel and endcap
                // electromagnetic calorimeters
                if (lep_type->at(i) == 11 && TMath::Abs(lep_eta->at(i)) < 2.47 &&
                    (TMath::Abs(lep_eta->at(i)) < 1.37 || TMath::Abs(lep_eta->at(i)) > 1.52)) {
                  if (TMath::Abs(lep_trackd0pvunbiased->at(i)) / lep_tracksigd0pvunbiased->at(i) < 5 &&
                      TMath::Abs(lep_z0->at(i) * TMath::Sin(leptemp.Theta())) < 0.5) {
                    goodlep_n     = goodlep_n + 1;
                    goodlep_index = i;
                    lep_index++;
                  }
                }
                // muon selection
                if (lep_type->at(i) == 13 && TMath::Abs(lep_eta->at(i)) < 2.5) {
                  if (TMath::Abs(lep_trackd0pvunbiased->at(i)) / lep_tracksigd0pvunbiased->at(i) < 3 &&
                      TMath::Abs(lep_z0->at(i) * TMath::Sin(leptemp.Theta())) < 0.5) {

                    goodlep_n     = goodlep_n + 1;
                    goodlep_index = i;
                    lep_index++;
                  }
                }
              }
            }
          }

          // Exactly one good lepton
          if (goodlep_n == 1) {

            // Preselection of good jets
            int goodjet_n  = 0;
            int goodbjet_n = 0;

            int goodjet_index[jet_n];
            int jet_index = 0;

            int goodbjet_index[jet_n];
            int bjet_index = 0;

            for (unsigned int i = 0; i < jet_n; i++) {
              if (jet_pt->at(i) > 30000. && TMath::Abs(jet_eta->at(i)) < 2.5) {
                // JVT cleaning
                bool jvt_pass = true;
                if (jet_pt->at(i) < 60000. && TMath::Abs(jet_eta->at(i)) < 2.4 && jet_jvt->at(i) < 0.59) jvt_pass = false;
                if (jvt_pass) {
                  goodjet_n++;
                  goodjet_index[jet_index] = i;
                  jet_index++;

                  // cut on 0.8244273 is 70% WP
                  if (jet_MV2c10->at(i) > 0.8244273) {
                    goodbjet_n++;
                    goodbjet_index[bjet_index] = i;
                    bjet_index++;
                  }
                }
              }
            }

            // TLorentzVector definitions
            TLorentzVector Lepton_1 = TLorentzVector();
            TLorentzVector MeT      = TLorentzVector();

            // nominal values
            Lepton_1.SetPtEtaPhiE(lep_pt->at(goodlep_index), lep_eta->at(goodlep_index), lep_phi->at(goodlep_index), lep_E->at(goodlep_index));
            MeT.SetPtEtaPhiE(met_et, 0, met_phi, met_et);

            // Calculation of MTW
            float mtw = sqrt(2.0f * Lepton_1.Pt() * MeT.Et() * (1.0f - cos(Lepton_1.DeltaPhi(MeT))));

            // At least four good jets
            if (goodjet_n >= 4) {
              int goodjet1_index = goodjet_index[0]; // leading jet

              // At least two b-tagged jets
              if (goodbjet_n == 2) {
                int goodbjet1_index = goodbjet_index[0]; // leading b-jet

                // MTW > 30 GeV
                if (mtw > 30000.) {
                  // systematic variations on objects

                  // First, we sample a random number from a Gaussian
                  // distribution with a given mean and sigma
                  TRandom3 *gRand            = new TRandom3(0);
                  Double_t  mean_met         = met_et;                           // nominal
                  Double_t  sigma_met        = met_et_syst;                      // width corresponds to the uncertainty
                  float     met_et_variation = gRand->Gaus(mean_met, sigma_met); // this is the variation up and down

                  // leptons syst
                  Double_t mean_lep         = lep_pt->at(goodlep_index);
                  Double_t sigma_lep        = lep_pt_syst->at(goodlep_index);
                  float    lep_pt_variation = gRand->Gaus(mean_lep, sigma_lep);

                  // lepton with syst
                  TLorentzVector Lepton_1_syst = TLorentzVector();
                  Lepton_1_syst.SetPtEtaPhiE(lep_pt_variation, lep_eta->at(goodlep_index), lep_phi->at(goodlep_index), lep_E->at(goodlep_index));

                  // MET with syst
                  TLorentzVector MeT_syst = TLorentzVector();
                  MeT_syst.SetPtEtaPhiE(met_et_variation, 0, met_phi, met_et_variation);

                  // MTW with syst variation on both MET and lepton pT
                  float mtw_syst = sqrt(2 * Lepton_1_syst.Pt() * MeT_syst.Et() * (1 - cos(Lepton_1_syst.DeltaPhi(MeT_syst))));

                  // syst variation on jet pT
                  Double_t mean_jet             = jet_pt->at(goodjet1_index) / 1000.;
                  Double_t sigma_jet            = (jet_pt_syst->at(goodjet1_index) / 1000.);
                  float    leadjet_pt_variation = gRand->Gaus(mean_jet, sigma_jet);

                  ///////// SAVE HISTOGRAMS /////////
                  double names_of_global_variable[] = {met_et / 1000., mtw / 1000.};

                  // correspond to the analysis cuts
                  if (MeT_syst.Et() / 1000. > 30 && mtw_syst / 1000. > 30) {
                    FillHistogramsGlobal(MeT_syst.Et() / 1000., weight, "hist_syst_etmiss");
                    FillHistogramsGlobal(mtw_syst / 1000., weight, "hist_syst_mtw");
                  }
                  if (Lepton_1_syst.Pt() / 1000. > 30) FillHistogramsLeadlept(Lepton_1_syst.Pt() / 1000., weight, "hist_syst_leadleptpt");

                  double names_of_leadlep_variable[] = {Lepton_1.Pt() / 1000.,
                                                        Lepton_1.Eta(),
                                                        Lepton_1.E() / 1000.,
                                                        Lepton_1.Phi(),
                                                        (double)lep_charge->at(goodlep_index),
                                                        (double)lep_type->at(goodlep_index),
                                                        Lepton_1.Et()/1000.0};

                  double names_of_jet_variable[] = {(double)goodjet_n,   jet_pt->at(goodjet1_index) / 1000.,  jet_eta->at(goodjet1_index),
                                                    (double)goodbjet_n,  jet_pt->at(goodbjet1_index) / 1000., jet_eta->at(goodbjet1_index),
                                                    leadjet_pt_variation};

                  TString histonames_of_global_variable[] = {"hist_etmiss", "hist_mtw"};

                  TString histonames_of_leadlep_variable[] = {"hist_leadleptpt", "hist_leadlepteta", "hist_leadleptE", "hist_leadleptphi",
                                                              "hist_leadleptch", "hist_leadleptID",  "hist_leadleptEt"};

                  TString histonames_of_jet_variable[] = {"hist_n_jets",      "hist_leadjet_pt",   "hist_leadjet_eta",    "hist_n_bjets",
                                                          "hist_leadbjet_pt", "hist_leadbjet_eta", "hist_syst_leadjet_pt"};

                  int length_global = sizeof(names_of_global_variable) / sizeof(names_of_global_variable[0]);

                  int length_leadlep = sizeof(names_of_leadlep_variable) / sizeof(names_of_leadlep_variable[0]);

                  int length_leadjet = sizeof(names_of_jet_variable) / sizeof(names_of_jet_variable[0]);

                  for (int i = 0; i < length_global; i++)
                    FillHistogramsGlobal(names_of_global_variable[i], weight, histonames_of_global_variable[i]);

                  for (int i = 0; i < length_leadlep; i++)
                    FillHistogramsLeadlept(names_of_leadlep_variable[i], weight, histonames_of_leadlep_variable[i]);

                  for (int i = 0; i < length_leadjet; i++)
                    FillHistogramsLeadJet(names_of_jet_variable[i], weight, histonames_of_jet_variable[i]);


                  // Invariant mass distribution of the 3-jets combination with
                  // the highest vector pT sum, a handle on the top mass

                  float PTjjjmax = 0;
                  float Mjjjmax  = 0;
                  int   a        = 0;
                  int   b        = 0;
                  int   c        = 0;

                  // iterate over 3 jets, build vectors
                  int bjet_jjj = -1;
                  //                  bool i_is_bjet, j_is_bjet, k_
                  for (int i = 0; i < goodjet_n; ++i) {
                    // int bjet_matched = 0;
                    int igji = goodjet_index[i];

                    for (int j = i + 1; j < goodjet_n; ++j) {
                      int igjj = goodjet_index[j];

                      for (int k = j + 1; k < goodjet_n; ++k) {
                        int igjk = goodjet_index[k];

                        int bjets_found = 0;
                        int bjet_used   = -1;
                        if (igji == goodbjet_index[0]) {
                          bjets_found++;
                          bjet_used = 0;
                        }
                        if (igjj == goodbjet_index[0]) {
                          bjets_found++;
                          bjet_used = 0;
                        }
                        if (igjk == goodbjet_index[0]) {
                          bjets_found++;
                          bjet_used = 0;
                        }
                        if (igji == goodbjet_index[1]) {
                          bjets_found++;
                          bjet_used = 1;
                        }
                        if (igjj == goodbjet_index[1]) {
                          bjets_found++;
                          bjet_used = 1;
                        }
                        if (igjk == goodbjet_index[1]) {
                          bjets_found++;
                          bjet_used = 1;
                        }
                        // if (!bjets_found) continue;

                        if (bjets_found != 1) continue;

                        TLorentzVector jet1 = TLorentzVector();
                        jet1.SetPtEtaPhiE(jet_pt->at(igji), jet_eta->at(igji), jet_phi->at(igji), jet_E->at(igji));
                        TLorentzVector jet2 = TLorentzVector();
                        jet2.SetPtEtaPhiE(jet_pt->at(igjj), jet_eta->at(igjj), jet_phi->at(igjj), jet_E->at(igjj));
                        TLorentzVector jet3 = TLorentzVector();
                        jet3.SetPtEtaPhiE(jet_pt->at(igjk), jet_eta->at(igjk), jet_phi->at(igjk), jet_E->at(igjk));

                        // find largest pT of 3-jet system,
                        float PTjjjTemp = (jet1 + jet2 + jet3).Pt() / 1000.;

                        if (PTjjjTemp > PTjjjmax) {
                          bjet_jjj = bjet_used; ///////

                          PTjjjmax = PTjjjTemp;
                          Mjjjmax  = (jet1 + jet2 + jet3).M() / 1000.; // this is m(jjj)
                          a        = i;
                          b        = j;
                          c        = k; // store the indices

                          // among those jets, find largest pT of 2-jet
                          // system, a handle of the W-boson
                          float PTjjTemp12 = (jet1 + jet2).Pt() / 1000.;
                          float PTjjTemp13 = (jet1 + jet3).Pt() / 1000.;
                          float PTjjTemp23 = (jet2 + jet3).Pt() / 1000.;

                          if (PTjjTemp12 > PTjjTemp13 && PTjjTemp12 > PTjjTemp23) {
                            a = i;
                            b = j;
                            c = k;
                          }
                          if (PTjjTemp13 > PTjjTemp12 && PTjjTemp13 > PTjjTemp23) {
                            a = i;
                            b = k;
                            c = j;
                          }
                          if (PTjjTemp23 > PTjjTemp12 && PTjjTemp23 > PTjjTemp13) {
                            a = j;
                            b = k;
                            c = i;
                          }
                        }
                      }
                    }
                  }

                  TLorentzVector bjet_1  = TLorentzVector();
                  int            bjet_lb = goodbjet_index[1 - bjet_jjj];
                  bjet_1.SetPtEtaPhiE(jet_pt->at(bjet_lb), jet_eta->at(bjet_lb), jet_phi->at(bjet_lb), jet_E->at(bjet_lb));

                  static const float m_t = 172.800f;
                  static const float M_W = 80.400f;
                  float Meb = (Lepton_1 + bjet_1).M() / 1000.0f;
                  costheta  = 2.0f * Meb * Meb / (m_t * m_t - M_W * M_W) - 1.0f;
                  
                  static const float zplus = -0.58740105f;
                  static const float zminus = 0.58740105f;


                  // FillHistogramsGlobal(costheta, 1.0, "hist_costheta_overflow");
                  hist_costheta_overflow->Fill(costheta, weight);
                  // Remove overflow & underflow
                  if (costheta >= -1.0f && costheta <= 1.0f) {
                    NTotal++;
                    if (costheta < zplus)  N1++;
                    if (costheta > zminus) N3++;

                    // FillHistogramsGlobal(costheta, 1.0, "hist_costheta");
                    hist_costheta->Fill(costheta, weight);
                  }


                  /////#define m_t 172800.0f
                  /////#define M_W 80400.0f
                  /////float MW2    = M_W*M_W;
                  /////float me     = Lepton_1.M();
                  /////float mb     = bjet_1.M();
                  /////float Ee     = (MW2 + me)/(2*M_W);
                  /////float Eb     = (m_t*m_t - MW2 - mb*mb)/(2*M_W);
                  /////float pe_mag = TMath::Sqrt(Ee*Ee - me*me);
                  /////float pb_mag = TMath::Sqrt(Eb*Eb - mb*mb);
                  /////costheta     = (Lepton_1.Dot(bjet_1) - Ee*Eb) / (pe_mag * pb_mag);

                  // among the previous 3 jets, we search for those 2 that have
                  // the maximum pT, a handle on the W-boson mass
                  TLorentzVector j1 = TLorentzVector();
                  j1.SetPtEtaPhiE(jet_pt->at(goodjet_index[a]), jet_eta->at(goodjet_index[a]), jet_phi->at(goodjet_index[a]),
                                  jet_E->at(goodjet_index[a]));
                  TLorentzVector j2 = TLorentzVector();
                  j2.SetPtEtaPhiE(jet_pt->at(goodjet_index[b]), jet_eta->at(goodjet_index[b]), jet_phi->at(goodjet_index[b]),
                                  jet_E->at(goodjet_index[b]));

                  float Mjjmax = (j1 + j2).M() / 1000.; // first indices

                  if (Mjjjmax > 100 && Mjjjmax < 250) FillHistogramsTTbar(Mjjjmax, weight, "hist_Topmass");

                  if (Mjjmax > 50 && Mjjmax < 120) FillHistogramsTTbar(Mjjmax, weight, "hist_Wmass");

                  // delete random number
                  gRand->Delete();

                  
                }
              }
            }
          }
        }
      }
    }
  }
  return kTRUE;
}

void TTbarAnalysis::SlaveTerminate() {}

void TTbarAnalysis::Terminate() {
  //int NTotal = hist_costheta->GetEntries();
  printf("NTotal : %d\n", NTotal);
  int N2 = NTotal - N1 - N3;
  double Aplus  = (double)(N2 +  N3 - N1)  / (double)(N1 + N2 + N3);
  double Aminus = (double)(N3 - (N2 + N1)) / (double)(N1 + N2 + N3);
  double beta = 0.25992104;
  double Fr = 1.0/(1.0-beta) + (Aminus - beta * Aplus)/(3.0*beta*(1.0-beta*beta));
  double Fl = 1.0/(1.0-beta) - (Aplus - beta * Aminus)/(3.0*beta*(1.0-beta*beta));
  double F0 = -(1.0+beta)/(1.0-beta) + (Aplus - Aminus)/(3.0*beta*(1.0-beta));

  printf("Value of Fr %g\n", Fr);
  printf("Value of Fl %g\n", Fl);
  printf("Value of F0 %g\n", F0);
  printf("Value of Aplus %g\n", Aplus);
  printf("Value of Aminus %g\n", Aminus);


  TString filename_option = GetOption();

  // root file name, Aplus, Aminus, Fr, Fl, F0, N1, N2, N3
  FILE* angular_asymmetries_file = fopen("angular_asymmetries_generated.csv", "a");
  fprintf(angular_asymmetries_file, "\n%s, \t%f, \t%f, ", (const char*)filename_option, Aplus, Aminus);
  fprintf(angular_asymmetries_file, "\t%f, \t%f, \t%f, ", Fr, Fl, F0);
  fprintf(angular_asymmetries_file, "\t%d, \t%d, \t%d", N1, N2, N3);
  fclose(angular_asymmetries_file);

  printf("Writing with name option: %s \n", filename_option.Data());
  TString     output_name = "Output_TTbarAnalysis/" + filename_option + ".root";
  const char *filename    = output_name;

  TFile physicsoutput_TTbar(filename, "recreate");
  WriteHistograms();
  physicsoutput_TTbar.Close();

  
}
