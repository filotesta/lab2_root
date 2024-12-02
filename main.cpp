#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TStyle.h"
#include "particle.hpp"
#include <TBenchmark.h>

Impulse polarToCartesian(double r, double theta, double phi)
{
  return {r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta)};
}

void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptTitle(0);
}

int generate()
{
  setStyle();
  TBenchmark benchmark;
  benchmark.Start("Execution");

  Particle::addParticleType("Pi+", 0.13957, 1, 0);
  Particle::addParticleType("Pi-", 0.13957, -1, 0);
  Particle::addParticleType("K+", 0.49367, 1, 0);
  Particle::addParticleType("K-", 0.49367, -1, 0);
  Particle::addParticleType("P+", 0.93827, 1, 0);
  Particle::addParticleType("P-", 0.93827, -1, 0);
  Particle::addParticleType("K*", 0.89166, 0, 0.05);

  gRandom->SetSeed();
  constexpr int nEvents{int(1e5)};
  constexpr int n{120};
  Particle eventParticles[n];
  int overflow{100};

  TString hTitle[12] = {"Particle Types",
                        "Angle Phi Distribution",
                        "Angle Theta Distribution",
                        "Impulse Distribution",
                        "Transverse Impulse Distribution",
                        "Particle Energy",
                        "Invariant Mass: all particles",
                        "Invariant Mass: particles with same charges",
                        "Invariant Mass: particles with opposite charges",
                        "Invariant Mass: Pi+/K+ or Pi-/K-",
                        "Invariant Mass: Pi-/K+ or Pi+/K-",
                        "Invariant mass: decayment particles"};
  TH1F* histo[12];
  histo[0] = new TH1F(TString("h0"), hTitle[0], 8, 0., 7.);
  histo[0]->GetYaxis()->SetTitleOffset(1.);
  histo[0]->GetXaxis()->SetTitleSize(0.05);
  histo[0]->GetXaxis()->CenterTitle(true);
  histo[0]->GetYaxis()->CenterTitle(true);
  histo[0]->GetYaxis()->SetTitleSize(0.05);
  histo[0]->GetXaxis()->SetTitle(hTitle[0]);
  histo[0]->GetYaxis()->SetTitle("Entries");
  histo[0]->SetLineColor(kBlue);
  histo[0]->SetMarkerColor(kBlue);
  histo[0]->SetMarkerStyle(7);
  histo[0]->SetLineWidth(1);

  for (int i{1}; i < 3; ++i) {
    histo[i] = new TH1F(TString("h") + i, hTitle[i], 200, 0., 2 * M_PI / i);
    histo[i]->GetYaxis()->SetTitleOffset(1.);
    histo[i]->GetXaxis()->SetTitleSize(0.05);
    histo[i]->GetXaxis()->CenterTitle(true);
    histo[i]->GetYaxis()->CenterTitle(true);
    histo[i]->GetYaxis()->SetTitleSize(0.05);
    histo[i]->GetXaxis()->SetTitle(hTitle[i]);
    histo[i]->GetYaxis()->SetTitle("Entries");
    histo[i]->SetLineColor(kBlue);
    histo[i]->SetMarkerColor(kBlue);
    histo[i]->SetMarkerStyle(7);
    histo[i]->SetLineWidth(1);
  }

  for (int i{3}; i < 6; ++i) {
    histo[i] = new TH1F(TString("h") + i, hTitle[i], 200, 0., 6.5);
    histo[i]->GetYaxis()->SetTitleOffset(1.);
    histo[i]->GetXaxis()->SetTitleSize(0.05);
    histo[i]->GetXaxis()->CenterTitle(true);
    histo[i]->GetYaxis()->CenterTitle(true);
    histo[i]->GetYaxis()->SetTitleSize(0.05);
    histo[i]->GetXaxis()->SetTitle(hTitle[i]);
    histo[i]->GetYaxis()->SetTitle("Entries");
    histo[i]->SetLineColor(kBlue);
    histo[i]->SetMarkerColor(kBlue);
    histo[i]->SetMarkerStyle(7);
    histo[i]->SetLineWidth(1);
  }

  for (int i{6}; i < 11; ++i) {
    histo[i] = new TH1F(TString("h") + i, hTitle[i], 3000, 0., 6.0);
    histo[i]->GetYaxis()->SetTitleOffset(1.);
    histo[i]->GetXaxis()->SetTitleSize(0.05);
    histo[i]->GetXaxis()->CenterTitle(true);
    histo[i]->GetYaxis()->CenterTitle(true);
    histo[i]->GetYaxis()->SetTitleSize(0.05);
    histo[i]->GetXaxis()->SetTitle(hTitle[i]);
    histo[i]->GetYaxis()->SetTitle("Entries");
    histo[i]->SetLineColor(kBlue);
    histo[i]->SetMarkerColor(kBlue);
    histo[i]->SetMarkerStyle(7);
    histo[i]->SetLineWidth(1);
  }

  histo[11] = new TH1F(TString("h11"), hTitle[11], 300, 0., 2.);
  histo[11]->GetYaxis()->SetTitleOffset(1.);
  histo[11]->GetXaxis()->SetTitleSize(0.05);
  histo[11]->GetXaxis()->CenterTitle(true);
  histo[11]->GetYaxis()->CenterTitle(true);
  histo[11]->GetYaxis()->SetTitleSize(0.05);
  histo[11]->GetXaxis()->SetTitle(hTitle[11]);
  histo[11]->GetYaxis()->SetTitle("Entries");
  histo[11]->SetLineColor(kBlue);
  histo[11]->SetMarkerColor(kBlue);
  histo[11]->SetMarkerStyle(7);
  histo[11]->SetLineWidth(1);

  for (int j{0}; j < nEvents; ++j) {
    overflow = 100;
    for (int i{0}; i < 100; ++i) {
      double phi     = gRandom->Uniform(0, 2 * M_PI);
      double theta   = gRandom->Uniform(0, M_PI);
      double impulse = gRandom->Exp(1);
      histo[1]->Fill(phi);
      histo[2]->Fill(theta);
      histo[3]->Fill(impulse);

      double x = gRandom->Uniform(0, 1);
      if (0. <= x && x < 0.4)
        eventParticles[i].setIndex("Pi+");
      else if (0.4 <= x && x < 0.8)
        eventParticles[i].setIndex("Pi-");
      else if (0.8 <= x && x < 0.85)
        eventParticles[i].setIndex("K+");
      else if (0.85 <= x && x < 0.9)
        eventParticles[i].setIndex("K-");
      else if (0.9 <= x && x < 0.945)
        eventParticles[i].setIndex("P+");
      else if (0.945 <= x && x < 0.99)
        eventParticles[i].setIndex("P-");
      else {
        eventParticles[i].setIndex("K*");
        Particle dau1{};
        Particle dau2{};
        if (x < 0.995) {
          dau1.setIndex("Pi+");
          dau2.setIndex("K-");
        } else {
          dau1.setIndex("Pi-");
          dau2.setIndex("K+");
        }

        if (eventParticles[i].Decay2Body(dau1, dau2) == 0) {
          eventParticles[overflow].setIndex(dau1.getIndex());
          ++overflow;
          eventParticles[overflow].setIndex(dau2.getIndex());
          ++overflow;
          histo[11]->Fill(dau1.particleInvMass(dau2));
        }
      }
      eventParticles[i].setImpulse(polarToCartesian(impulse, theta, phi));
      histo[4]->Fill(
          sqrt(std::pow(impulse, 2) + std::pow(eventParticles[i].getImpulse().px_, 2)));
      histo[0]->Fill(eventParticles[i].getIndex());
      histo[5]->Fill(eventParticles[i].particleEnergy());
    }

    for (int k{0}; k < overflow - 1; ++k) {
      for (int r{k + 1}; r < overflow; ++r) {
        if (eventParticles[k].getIndex() != 6 && eventParticles[r].getIndex() != 6) {
          histo[6]->Fill(eventParticles[k].particleInvMass(eventParticles[r]));
        }
      }
    } // massa invariante fra tutte le particelle

    for (int k{0}; k < overflow - 1; ++k) {
      for (int r{k + 1}; r < overflow; ++r) {
        if (eventParticles[k].getIndex() != 6 && eventParticles[r].getIndex() != 6) {
          if (eventParticles[k].getCharge() * eventParticles[r].getCharge() == 1) { // carica concorde
            histo[7]->Fill(eventParticles[k].particleInvMass(eventParticles[r]));
          } else {
            assert((eventParticles[k].getCharge() * eventParticles[r].getCharge()) == -1); // carica discorde
            histo[8]->Fill(eventParticles[k].particleInvMass(eventParticles[r]));
          }
        }
      }
    }

    for (int k{0}; k < overflow - 1; ++k) {
      for (int r{k + 1}; r < overflow; ++r) {
        if ((eventParticles[k].getIndex() + eventParticles[r].getIndex()) == 3) {
          histo[10]->Fill(eventParticles[k].particleInvMass(eventParticles[r]));
        } else if ((eventParticles[k].getIndex() < 4 && eventParticles[r].getIndex() < 4)
                   && eventParticles[k].getIndex() != eventParticles[r].getIndex()
                   && ((eventParticles[k].getIndex() + eventParticles[r].getIndex() == 2)
                       || (eventParticles[k].getIndex() + eventParticles[r].getIndex()
                           == 4))) {
          histo[9]->Fill(eventParticles[k].particleInvMass(eventParticles[r]));
        }
      }
    }
  }

  for (int i{7}; i < 11; ++i) {
    histo[i]->Sumw2();
  }

  TCanvas* canvas1 = new TCanvas(
      "canvas1", "Informazioni particelle: tipo, angoli, impulsi ed energia", 0, 10, 800, 600);
  canvas1->Divide(2, 3);
  for (int i{0}; i < 6; ++i) {
    canvas1->cd(i + 1);
    histo[i]->Draw("E");
    histo[i]->Draw("HISTO,SAME");
  }

  TCanvas* canvas2 = new TCanvas("canvas2", "Informazioni particelle: massa invariante", 1000,
                                 1000, 800, 600);
  canvas2->Divide(2, 3);
  for (int i{6}, j{0}; i < 12; ++i, ++j) {
    canvas2->cd(j + 1);
    histo[i]->Draw("E");
    histo[i]->Draw("HISTO,SAME");
  }

  benchmark.Stop("Execution");
  benchmark.Show("Execution");

  TFile* file = new TFile("particle_simulation.root", "RECREATE");
  for (int i{0}; i < 12; ++i) {
    histo[i]->Write();
  }
  file->Close();
  return 0;
}