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

  //Generazone delle canvas che verranno utilizzate
  TCanvas* canvasAngImp = new TCanvas(
      "canvasAngImp", "Fitting histograms about angles and impulses", 0, 10, 800, 600);

  TCanvas* canvasDecay =
      new TCanvas("canvasDecay", "Fitting histograms about decayments", 0, 10, 800, 600);

  // Fit: angoli phi, theta e impulso
  gStyle->SetOptFit(1111);
  TF1* phiUniformFit = new TF1("phiUniformFit", "[0]", hPhi->GetXaxis()->GetXmin(),
                               hPhi->GetXaxis()->GetXmax());
  phiUniformFit->SetLineColor(kRed);
  std::cout << "\n------------------------------------------------------------------------\n";
  hPhi->GetXaxis()->SetTitle("");
  hPhi->Fit("phiUniformFit", "R");
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

  TF1* thetaUniformFit = new TF1("thetaUniformFit", "[0]", hTheta->GetXaxis()->GetXmin(),
                                 hTheta->GetXaxis()->GetXmax());
  thetaUniformFit->SetLineColor(kRed);
  hTheta->GetXaxis()->SetTitle("");
  hTheta->Fit("thetaUniformFit", "R");
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

  TF1* impulseExpFit = new TF1("impulseExpFit", "expo", hImpulse->GetXaxis()->GetXmin(),
                               hImpulse->GetXaxis()->GetXmax());
  impulseExpFit->SetLineColor(kRed);
  hImpulse->GetXaxis()->SetTitle("");
  hImpulse->Fit("impulseExpFit", "R");
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

  // Generazione istogrammi per le particelle + fit
  TH1F* hSubtraction1 = new TH1F(
      "hS1", "Invariant Mass: subtraction between opposite and same charge", 1200, 0., 6.);
  hSubtraction1->Add(hInvMassOppositeCharge, hInvMassSameCharge, 1., -1.);
  hSubtraction1->GetYaxis()->SetTitleOffset(1.);
  hSubtraction1->GetXaxis()->SetTitleSize(0.05);
  hSubtraction1->GetXaxis()->CenterTitle(true);
  hSubtraction1->GetYaxis()->CenterTitle(true);
  hSubtraction1->GetYaxis()->SetTitleSize(0.05);
  hSubtraction1->GetXaxis()->SetTitle("");
  hSubtraction1->GetYaxis()->SetTitle("Entries");
  hSubtraction1->SetLineColor(kBlue);
  hSubtraction1->SetMarkerColor(kBlue);
  hSubtraction1->SetMarkerStyle(7);
  hSubtraction1->SetLineWidth(1);

  TH1F* hSubtraction2 =
      new TH1F("hS2", "Invariant Mass: subtraction between Pi and K", 1200, 0., 6.);
  hSubtraction2->Add(hInvMassPiKOppositeCharge, hInvMassPiKSameCharge, 1., -1.);
  hSubtraction2->GetYaxis()->SetTitleOffset(1.);
  hSubtraction2->GetXaxis()->SetTitleSize(0.05);
  hSubtraction2->GetXaxis()->CenterTitle(true);
  hSubtraction2->GetYaxis()->CenterTitle(true);
  hSubtraction2->GetYaxis()->SetTitleSize(0.05);
  hSubtraction2->GetXaxis()->SetTitle("");
  hSubtraction2->GetYaxis()->SetTitle("Entries");
  hSubtraction2->SetLineColor(kBlue);
  hSubtraction2->SetMarkerColor(kBlue);
  hSubtraction2->SetMarkerStyle(7);
  hSubtraction2->SetLineWidth(1);

  TF1* fitGauss1 = new TF1("fitGauss1", "gaus", hSubtraction1->GetXaxis()->GetXmin(),
                           hSubtraction1->GetXaxis()->GetXmax());
  fitGauss1->SetLineColor(kRed);
  fitGauss1->SetParameter(1, 0.89166);
  fitGauss1->SetParameter(2, 0.05);
  hSubtraction1->Fit(fitGauss1, "R");
  hSubtraction1->GetXaxis()->SetRangeUser(0.2, 1.5);
  std::cout << "\n                        Fit hSubtraction1                        \n"
            << "Mean: " << fitGauss1->GetParameter(1) << "\n";
  std::cout << "Sigma: " << fitGauss1->GetParameter(2) << "\n";
  std::cout << "Recudced Chisquare: " << fitGauss1->GetChisquare() / fitGauss1->GetNDF()
            << "\n";
  std::cout << "Fit probability = " << fitGauss1->GetProb() << "\n";
  std::cout << "\n------------------------------------------------------------------------\n";

  TF1* fitGauss2 = new TF1("fitGauss2", "gaus", hSubtraction2->GetXaxis()->GetXmin(),
                           hSubtraction2->GetXaxis()->GetXmax());
  fitGauss2->SetLineColor(kRed);
  fitGauss2->SetParameter(1, 0.89166);
  fitGauss2->SetParameter(2, 0.05);
  hSubtraction2->Fit(fitGauss2, "R");
  hSubtraction2->GetXaxis()->SetRangeUser(0.2, 1.5);
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

  canvasAngImp->Divide(1, 3);
  canvasAngImp->cd(1);
  hPhi->Draw();
  canvasAngImp->cd(2);
  hTheta->Draw();
  canvasAngImp->cd(3);
  hImpulse->Draw();

  canvasDecay->Divide(1, 3);
  canvasDecay->cd(1);
  hInvMassKStar->GetXaxis()->SetTitle("");
  hInvMassKStar->Draw();
  canvasDecay->cd(2);
  hSubtraction1->Draw();
  canvasDecay->cd(3);
  hSubtraction2->Draw();

  TFile* file2 = new TFile("particle_checking.root", "RECREATE");
  hInvMassKStar->Write();
  hSubtraction2->Write();
  hSubtraction1->Write();
  file2->Close();

  canvasAngImp->Print("canvasAngImp.pdf", "RECREATE");
  canvasDecay->Print("canvasDecay.pdf", "RECREATE");
}
