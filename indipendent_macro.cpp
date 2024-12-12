#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TStyle.h"
#include "particle.hpp"
#include <iostream>

void analysis()
{
  TFile* file                     = TFile::Open("particle_simulation.root", "READ");
  TH1F* hParticles                = (TH1F*)file->Get(TString("h0"));
  TH1F* hPhi                      = (TH1F*)file->Get(TString("h1"));
  TH1F* hTheta                    = (TH1F*)file->Get(TString("h2"));
  TH1F* hImpulse                  = (TH1F*)file->Get(TString("h3"));
  TH1F* hTransvImpulse            = (TH1F*)file->Get(TString("h4"));
  TH1F* hEnergy                   = (TH1F*)file->Get(TString("h5"));
  TH1F* hInvMass                  = (TH1F*)file->Get(TString("h6"));
  TH1F* hInvMassSameCharge        = (TH1F*)file->Get(TString("h7"));
  TH1F* hInvMassOppositeCharge    = (TH1F*)file->Get(TString("h8"));
  TH1F* hInvMassPiKSameCharge     = (TH1F*)file->Get(TString("h9"));
  TH1F* hInvMassPiKOppositeCharge = (TH1F*)file->Get(TString("h10"));
  TH1F* hInvMassKStar             = (TH1F*)file->Get(TString("h11"));

  const double nEntrieshParticle = hParticles->GetEntries();
  const double nEntrieshPhi      = hPhi->GetEntries();
  const double nEntrieshTheta    = hTheta->GetEntries();

  TCanvas* canvasPhi = new TCanvas("canvasPhi", "Fitting angle phi", 1000, 1000, 800, 600);

  gStyle->SetOptFit(1111);

  TF1* phiUniformFit = new TF1("phiUniformFit", "[0]", hPhi->GetXaxis()->GetXmin(),
                               hPhi->GetXaxis()->GetXmax());
  std::cout << "\n------------------------------------------------------------------------\n";
  hPhi->Fit("phiUniformFit");
  hPhi->Draw();
  std::cout << "\n                        Fit hPhi                        \n> Parameter: "
            << phiUniformFit->GetParameter(0) << "\n> Expected Value: "
            << nEntrieshPhi * (hPhi->GetXaxis()->GetXmax() - hPhi->GetXaxis()->GetXmin())
                   / (hPhi->GetNbinsX() * ((2 * M_PI)))
            << "\n> Parameter - Expected value: "
            << std::abs(phiUniformFit->GetParameter(0)
                        - nEntrieshPhi
                              * (hPhi->GetXaxis()->GetXmax() - hPhi->GetXaxis()->GetXmin())
                              / (hPhi->GetNbinsX() * ((2 * M_PI))))
            << "\n\n";
  std::cout << "> Reduced Chisquare: "
            << phiUniformFit->GetChisquare() / phiUniformFit->GetNDF();
  std::cout << "\n> Fit Probability: " << phiUniformFit->GetProb() << "\n";
  std::cout << "\n------------------------------------------------------------------------\n";

  TCanvas* canvasTheta =
      new TCanvas("canvasTheta", "Fitting angle theta", 1000, 1000, 800, 600);
  TF1* thetaUniformFit = new TF1("thetaUniformFit", "[0]", hTheta->GetXaxis()->GetXmin(),
                                 hTheta->GetXaxis()->GetXmax());
  hTheta->Fit("thetaUniformFit", "R");
  hTheta->Draw();
  std::cout << "\n                        Fit hTheta                        \n> Parameter: "
            << thetaUniformFit->GetParameter(0) << "\n> Expected Value: "
            << nEntrieshTheta * (hTheta->GetXaxis()->GetXmax() - hTheta->GetXaxis()->GetXmin())
                   / (hTheta->GetNbinsX() * ((M_PI)))
            << "\n> Parameter - Expected value: "
            << std::abs(thetaUniformFit->GetParameter(0)
                        - nEntrieshTheta
                              * (hTheta->GetXaxis()->GetXmax() - hTheta->GetXaxis()->GetXmin())
                              / (hTheta->GetNbinsX() * ((M_PI))))
            << "\n";
  std::cout << "> Reduced Chisquare: "
            << thetaUniformFit->GetChisquare() / thetaUniformFit->GetNDF();
  std::cout << "\n> Fit Probability: " << thetaUniformFit->GetProb() << "\n";
  std::cout << "\n------------------------------------------------------------------------\n";

  TCanvas* canvasImpulse = new TCanvas("canvasImpulse", "Fitting impulse", 1000, 1000, 800, 600);
  TF1* impulseExpFit = new TF1("impulseExpFit", "expo", hImpulse->GetXaxis()->GetXmin(),
                               hImpulse->GetXaxis()->GetXmax());
  hImpulse->Fit("impulseExpFit", "R");
  hImpulse->Draw();
  std::cout << "\n                        Fit hImpulse                        \n"
            << "> Mean: " << hImpulse->GetMean() << " +/- " << hImpulse->GetMeanError()
            << "\n";
  const double meanDiff = std::abs(hImpulse->GetMean() - 1);
  if (meanDiff <= hImpulse->GetMeanError()) {
    std::cout << "> Mean is COMPATIBLE with 1, with a difference of: " << meanDiff << "\n";
  } else {
    std::cout << "> Mean is NOT COMPATIBLE with 1, with a difference of: " << meanDiff << "\n";
  }

  std::cout << "> Parameter constant: " << impulseExpFit->GetParameter(0) << "\n"
            << "> Parameter slope: " << impulseExpFit->GetParameter(1) << "\n";
  std::cout << "> Reduced Chisquare: "
            << impulseExpFit->GetChisquare() / impulseExpFit->GetNDF() << "\n";
  std::cout << "> Fit probability: " << impulseExpFit->GetProb() << "\n";
  std::cout << "\n------------------------------------------------------------------------\n";

  TCanvas* canvasKStar = new TCanvas("canvasDecayment", "Fitting histogram from decayment K*",
                                     1000, 1000, 800, 600);
  hInvMassKStar->Draw();

  TH1F* hSubtraction1 = new TH1F(
      "hS1", "Invariant Mass: subtraction between opposite and same charge", 1200, 0., 6.);
  hSubtraction1->Add(hInvMassOppositeCharge, hInvMassSameCharge, 1., -1.);

  TH1F* hSubtraction2 =
      new TH1F("hS2", "Invariant Mass: subtraction between Pi and K", 1200, 0., 6.);
  hSubtraction2->Add(hInvMassPiKOppositeCharge, hInvMassPiKSameCharge, 1., -1.);

  TCanvas* canvasSub1 = new TCanvas(
      "canvasSub1", "Fitting histogram subtraction between opposite and same charge", 1000,
      1000, 800, 600);

  TF1* fitGauss1 = new TF1("fitGauss1", "gaus", hSubtraction1->GetXaxis()->GetXmin(),
                           hSubtraction1->GetXaxis()->GetXmax());
  fitGauss1->SetParameter(1, hInvMassKStar->GetMean());
  fitGauss1->SetParameter(2, hInvMassKStar->GetRMS());
  hSubtraction1->Fit(fitGauss1, "R");
  hSubtraction1->GetXaxis()->SetRangeUser(0.2, 1.5);
  hSubtraction1->Draw();
  std::cout << "\n                        Fit hSubtraction1                        \n"
            << "Mean: " << fitGauss1->GetParameter(1) << "\n";
  std::cout << "Sigma: " << fitGauss1->GetParameter(2) << "\n";
  std::cout << "Recudced Chisquare: " << fitGauss1->GetChisquare() / fitGauss1->GetNDF()
            << "\n";
  std::cout << "Fit probability = " << fitGauss1->GetProb() << "\n";
  std::cout << "\n------------------------------------------------------------------------\n";

 TCanvas* canvasSub2 = new TCanvas(
      "canvasSub2", "Fitting histogram subtraction between Pi and K", 1000,
      1000, 800, 600);
  TF1* fitGauss2 = new TF1("fitGauss2", "gaus", hSubtraction2->GetXaxis()->GetXmin(),
                           hSubtraction2->GetXaxis()->GetXmax());
  fitGauss2->SetParameter(1, hInvMassKStar->GetMean());
  fitGauss2->SetParameter(2, hInvMassKStar->GetRMS());
  hSubtraction2->Fit(fitGauss2, "R");
  hSubtraction2->GetXaxis()->SetRangeUser(0.2, 1.5);
  hSubtraction2->Draw();
  std::cout << "\n                        Fit hSubtraction2                        \n"
            << "Mean: " << fitGauss2->GetParameter(1) << "\n";
  std::cout << "Sigma: " << fitGauss2->GetParameter(2) << "\n";
  std::cout << "Recudced Chisquare: " << fitGauss2->GetChisquare() / fitGauss2->GetNDF()
            << "\n";
  std::cout << "Fit probability = " << fitGauss2->GetProb() << "\n";
  std::cout << "\n------------------------------------------------------------------------\n";

  std::array<std::string, 7> particleNames{"Pi+", "Pi-", "K+", "K-", "P+", "P-", "K*"};
  std::cout << "\n                        Particle Generation Ratios                        "
            << "\n> Total entries: " << nEntrieshParticle << "\n";
  for (int i{}; i < 7; ++i) {
    std::cout << "> " << particleNames[i] << " entries: " << hParticles->GetBinContent(i + 1)
              << " +/- " << hParticles->GetBinError(i + 1) << " ("
              << hParticles->GetBinContent(i + 1) * 100. / nEntrieshParticle << "%)\n";
  }

  TFile* file2 = new TFile("particle_checking.root", "RECREATE");
  hInvMassKStar->Write();
  hSubtraction2->Write();
  hSubtraction1->Write();
  file2->Close();

  canvasPhi->Print("canvasPhi.pdf","RECREATE");
  canvasTheta->Print("canvasTheta.pdf","RECREATE");
  canvasImpulse->Print("canvasImpulse.pdf","RECREATE");
  canvasKStar->Print("canvasKStar.pdf","RECREATE");
  canvasSub1->Print("canvasSub1.pdf","RECREATE");
  canvasSub2->Print("canvasSub2.pdf","RECREATE");

}
