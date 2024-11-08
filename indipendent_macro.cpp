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

  const double nEntrieshParticle       = hParticles->GetEntries();
  const double nEntrieshPhi            = hPhi->GetEntries();
  const double nEntrieshTheta          = hTheta->GetEntries();
  const double nEntrieshImpulse        = hImpulse->GetEntries();
  const double nEntrieshTransvImpulsee = hTransvImpulse->GetEntries();
  const double nEntrieshEnergy         = hEnergy->GetEntries();
  const double nEntrieshInvMass        = hInvMass->GetEntries();

  for (int bin = 1; bin < hParticles->GetNbinsX(); ++bin) {
    double content = hParticles->GetBinContent(bin);
    double error   = hParticles->GetBinError(bin);
  }

  TF1* phiUniformFit = new TF1("phiUniformFit", "[0]", hPhi->GetXaxis()->GetXmin(),
                               hPhi->GetXaxis()->GetXmax());
  hPhi->Fit("phiUniformFit");
  std::cout << "\n> Fit hPhi:\n> Parameter: " << phiUniformFit->GetParameter(0)
            << "\n> Expected Value: "
            << nEntrieshPhi * (hPhi->GetXaxis()->GetXmax() - hPhi->GetXaxis()->GetXmin())
                   / (hPhi->GetNbinsX() * ((2 * M_PI)))
            << "\n> Parameter - Expected value: "
            << std::abs(phiUniformFit->GetParameter(0)
                        - nEntrieshPhi
                              * (hPhi->GetXaxis()->GetXmax() - hPhi->GetXaxis()->GetXmin())
                              / (hPhi->GetNbinsX() * ((2 * M_PI))))
            << "\n\n";
  // std::cout << "> Reduced Chisquare: "
  //           << phiUniformFit->GetChisquare() / phiUniformFit->GetNDF();
  // std::cout << "\n> Fit Probability: " << phiUniformFit->GetProb() << "\n";

  TF1* thetaUniformFit = new TF1("thetaUniformFit", "[0]", hTheta->GetXaxis()->GetXmin(),
                               hTheta->GetXaxis()->GetXmax());
  hTheta->Fit("thetaUniformFit");
  std::cout << "\n> Fit hTheta:\n> Parameter: " << thetaUniformFit->GetParameter(0)
            << "\n> Expected Value: "
            << nEntrieshTheta * (hTheta->GetXaxis()->GetXmax() - hTheta->GetXaxis()->GetXmin())
                   / (hTheta->GetNbinsX() * ((M_PI)))
            << "\n> Parameter - Expected value: "
            << std::abs(thetaUniformFit->GetParameter(0)
                        - nEntrieshTheta
                              * (hTheta->GetXaxis()->GetXmax() - hTheta->GetXaxis()->GetXmin())
                              / (hTheta->GetNbinsX() * ((M_PI))))
            << "\n\n";
  // std::cout << "> Reduced Chisquare: "
  //           << thetaUniformFit->GetChisquare() / thetaUniformFit->GetNDF();
  // std::cout << "\n> Fit Probability: " << thetaUniformFit->GetProb() << "\n";

  TF1* impulseExpFit =
      new TF1("impulseExpFit", "expo", hImpulse->GetXaxis()->GetXmin(),
              hImpulse->GetXaxis()->GetXmax());
  hImpulse->Fit("impulseExpFit");
  std::cout << "\n> Fit hImpulse: \n"
            << "> Mean: " << hImpulse->GetMean() << " +/- " << hImpulse->GetMeanError() << "\n";
  const double meanDiff = std::abs(hImpulse->GetMean() - 1);
  if (meanDiff <= hImpulse->GetMeanError()) {
    std::cout << "> Mean is COMPATIBLE with 1, with a difference of: " << meanDiff << "\n";
  } else {
    std::cout << "> Mean is NOT COMPATIBLE with 1, with a difference of: " << meanDiff << "\n";
  }

  std::cout << "> Parameter constant: " << impulseExpFit->GetParameter(0) << "\n"
            << "> Parameter slope: " << impulseExpFit->GetParameter(1) << "\n\n";
  // std::cout << "Reduced Chisquare: "
  //           << impulseExpFit->GetChisquare() / impulseExpFit->GetNDF();
  // std::cout << "Fit probability = " << impulseExpFit->GetProb() << "\n";

  TH1F* hSubtraction1 = new TH1F(
      "hS1", "Invariant Mass: subtraction between opposite and same charge", 1000, -0.3, 6.5);
  hSubtraction1->Add(hInvMassOppositeCharge, hInvMassSameCharge, 1, -1);

  TH1F* hSubtraction2 =
      new TH1F("hS2", "Invariant Mass: subtraction between Pi and K", 1000, -0.3, 6.5);
  hSubtraction2->Add(hInvMassPiKOppositeCharge, hInvMassPiKSameCharge, 1, -1);

  // disegno per comparare i picchi
  TCanvas* pickComparison = new TCanvas("pickComparison", "Invariant Mass Analisys", 800, 600);
  pickComparison->Divide(1, 3);
  pickComparison->cd(1);
  hInvMassKStar->Draw();

  TF1* fitGauss1 = new TF1("fitGauss1", "gaus", hSubtraction1->GetXaxis()->GetXmin(),
                           hSubtraction1->GetXaxis()->GetXmax());

  fitGauss1->SetParameter(1, hInvMassKStar->GetMean());
  fitGauss1->SetParameter(2, hInvMassKStar->GetRMS());
  hSubtraction1->Fit("fitGauss1");
  pickComparison->cd(2);
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


  TFile* file2 = new TFile("particle_checking.root", "RECREATE");
  hInvMassKStar->Write();
  hSubtraction2->Write();
  hSubtraction1->Write();
  file2->Close();

  pickComparison->Print("CheckHistogram.pdf");
  // hSubtraction2->Draw("SAME");
}
