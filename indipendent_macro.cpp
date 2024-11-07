#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <iostream>

void CheckHistogram()
{
  TFile* file = TFile::Open("particle_simulation.root", "READ");

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

  double nEntrieshParticle       = hParticles->GetEntries();
  double nEntrieshPhi            = hPhi->GetEntries();
  double nEntrieshTheta          = hTheta->GetEntries();
  double nEntrieshImpulse        = hImpulse->GetEntries();
  double nEntrieshTransvImpulsee = hTransvImpulse->GetEntries();
  double nEntrieshEnergy         = hEnergy->GetEntries();
  double nEntrieshInvMass        = hInvMass->GetEntries();

  for (int bin = 1; bin < hParticles->GetNbinsX(); ++bin) {
    double content = hParticles->GetBinContent(bin);
    double error   = hParticles->GetBinError(bin);
  }

  TF1* phiUniformFit = new TF1("phiUniformFit", "1/([1]-[0])", hPhi->GetXaxis()->GetXmin(),
                               hPhi->GetXaxis()->GetXmax());
  hPhi->Fit("phiUniformFit");
  std::cout << "\n Fit hPhi:\n " << "Parameter a: " << phiUniformFit->GetParameter(0) << "\n"
            << "Parameter b: " << phiUniformFit->GetParameter(1) << "\n";
  std::cout << "Parameter a  fit - expected: " << phiUniformFit->GetParameter(0) << "\n"
            << "Parameter b fit - expected: " << phiUniformFit->GetParameter(1) - 2 * M_PI
            << "\n";
  std::cout << "Reduced Chisquare: "
            << phiUniformFit->GetChisquare() / phiUniformFit->GetNDF();
  std::cout << "Probabilità del fit = " << phiUniformFit->GetProb() << "\n";

  TF1* thetaUniformFit = new TF1("thetaUniformFit", "1/([1]-[0])",
                                 hTheta->GetXaxis()->GetXmin(), hTheta->GetXaxis()->GetXmax());
  hTheta->Fit("thetaUniformFit");
  std::cout << "\n \n Fit hPhi:\n " << "Parameter a: " << thetaUniformFit->GetParameter(0)
            << "\n"
            << "Parameter b: " << thetaUniformFit->GetParameter(1) << "\n";
  std::cout << "Parameter a  fit - expected: " << thetaUniformFit->GetParameter(0) << "\n"
            << "Parameter b fit - expected: " << thetaUniformFit->GetParameter(1) - M_PI
            << "\n";
  std::cout << "Reduced Chisquare: "
            << thetaUniformFit->GetChisquare() / thetaUniformFit->GetNDF();
  std::cout << "Fit probability = " << thetaUniformFit->GetProb() << "\n";

  TF1* impulseExpFit =
      new TF1("impulseExpFit", "[0]*exp(-x/[1])", hImpulse->GetXaxis()->GetXmin(),
              hImpulse->GetXaxis()->GetXmax());
  hImpulse->Fit("impulseExpFit");
  std::cout << "\n\n Fit hImpulse: \n"
            << "Mean: " << hImpulse->GetMean() << " +/- " << hImpulse->GetMeanError() << "\n";
  double meanDiff = hImpulse->GetMean() - 1;
  if (meanDiff <= hImpulse->GetMeanError()) {
    std::cout << "Mean is COMPATIBLE with 1, with a difference of: " << meanDiff << "\n";
  } else {
    std::cout << "Mean is INCOMPATIBLE with 1, with a difference of: " << meanDiff << "\n";
  }

  std::cout << "Parameter c: " << impulseExpFit->GetParameter(0) << "\n"
            << "Parameter tau: " << impulseExpFit->GetParameter(1) << "\n";
  std::cout << "Reduced Chisquare: "
            << impulseExpFit->GetChisquare() / impulseExpFit->GetNDF();
  std::cout << "Fit probability = " << impulseExpFit->GetProb() << "\n";

  hInvMassSameCharge->Sumw2();
  hInvMassOppositeCharge->Sumw2();
  hInvMassPiKSameCharge->Sumw2();
  hInvMassPiKOppositeCharge->Sumw2();

  TH1F* hSubtraction1 = new TH1F(
      "hS1", "Invariant Mass: subtraction between opposite and same charge", 1000, -0.3, 6.5);
  hSubtraction1->Add(hInvMassOppositeCharge, hInvMassSameCharge, 1, -1);

  TH1F* hSubtraction2 =
      new TH1F("hS2", "Invariant Mass: subtraction between Pi and K", 1000, -0.3, 6.5);
  hSubtraction2->Add(hInvMassPiKOppositeCharge, hInvMassPiKSameCharge, 1, -1);

  // disegno per comparare i picchi
  TCanvas* pickComparison = new TCanvas("pickComparison", "Invariant Mass Analisys", 800, 600);
  pickComparison->Divide(1, 3);
  pickComparison->cd(2);
  hInvMassKStar->Draw();

  TF1* fitGauss1 = new TF1("fitGauss1", "gaus", hSubtraction1->GetXaxis()->GetXmin(),
                           hSubtraction1->GetXaxis()->GetXmax());

  fitGauss1->SetParameter(1, hInvMassKStar->GetMean());
  fitGauss1->SetParameter(2, hInvMassKStar->GetRMS());
  hSubtraction1->Fit("fitGauss1");
  pickComparison->cd(1);
  hSubtraction1->Draw();
  std::cout << "\n Fit hSubtraction1: \n" << "Mean: " << fitGauss1->GetParameter(1) << "\n";
  std::cout << "Sigma: " << fitGauss1->GetParameter(2) << "\n";
  std::cout << "Recudced Chisquare: " << fitGauss1->GetChisquare() / fitGauss1->GetNDF()
            << "\n";
  std::cout << "Fit probability = " << fitGauss1->GetProb() << "\n";

  // hSubtraction1->Draw("SAME");

  TF1* fitGauss2 = new TF1("fitGauss2", "gaus", hSubtraction2->GetXaxis()->GetXmin(),
                           hSubtraction2->GetXaxis()->GetXmax());

  fitGauss2->SetParameter(1, hInvMassKStar->GetMean());
  fitGauss2->SetParameter(2, hInvMassKStar->GetRMS());
  hSubtraction2->Fit("fitGauss2");
  pickComparison->cd(3);
  hSubtraction2->Draw();
  std::cout << "\n Fit hSubtraction2: \n" << "Mean: " << fitGauss2->GetParameter(1) << "\n";
  std::cout << "Sigma: " << fitGauss2->GetParameter(2) << "\n";
  std::cout << "Recudced Chisquare: " << fitGauss2->GetChisquare() / fitGauss2->GetNDF()
            << "\n";
  std::cout << "Fit probability = " << fitGauss2->GetProb() << "\n";

  // hSubtraction2->Draw("SAME");
}
