R__LOAD_LIBRARY("/lab2_root/particle_cpp.so");
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TStyle.h"
#include "particle.hpp"

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

int main()
{
  setStyle();
  Particle::addParticleType("Pi+", 0.13957, 1, 0);
  Particle::addParticleType("Pi-", 0.13957, -1, 0);
  Particle::addParticleType("K+", 0.49367, 1, 0);
  Particle::addParticleType("K-", 0.49367, -1, 0);
  Particle::addParticleType("P+", 0.93827, 1, 0);
  Particle::addParticleType("P-", 0.93827, -1, 0);
  Particle::addParticleType("K*", 0.89166, 0, 0.05);

  gRandom->SetSeed();
  constexpr int nEvents{int(10e5)};
  constexpr int n{120};
  Particle eventParticles[n];
  int overflow{100};

  TH1F* hParticles     = new TH1F("hParticles", "particle types", 9, -1., 8.);
  TH1F* hPhi           = new TH1F("hPhi", "angle phi", 1000, -1, 2 * M_PI + 1);
  TH1F* hTheta         = new TH1F("hTheta", "angle theta", 1000, -1, M_PI + 1);
  TH1F* hImpulse       = new TH1F("hImpulse", "impulse", 1000, -1, 6);
  TH1F* hTransvImpulse = new TH1F("hTrasvImpulse", "trasverse impulse", 1000, -1, 6);
  TH1F* hEnergy        = new TH1F("hEnergy", "particle energy", 1000, -1, 6);
  TH1F* hInvMass       = new TH1F("hInvMass", "invariant mass", 1000, -1, 3);
  TH1F* hInvMass1 = new TH1F("hInvMass1", "invariant mass between particles with charges of the same sign", 1000, -1, 3);
  TH1F* hInvMass2 = new TH1F("hInvMass2", "invariant mass between particles with charges of different sign", 1000, -1, 3);
  TH1F* hInvMass3 = new TH1F("hInvMass3", "invariant mass between Pi and K with charges of different sign", 1000, -1, 3);
  TH1F* hInvMass4 = new TH1F("hInvMass4", "invariant mass between Pi and K with charges of the same sign", 1000, -1, 3);
  TH1F* hInvMass5 = new TH1F("hInvMass5", "invariant mass between particles decadeted from the same K*", 1000,-1, 3);

    for (int j{0}; j < nEvents; ++j) {
    for (int i{0}; i < 100; ++i) {
      double phi     = gRandom->Uniform(0, 2 * M_PI);
      double theta   = gRandom->Uniform(0, M_PI);
      double impulse = gRandom->Exp(-1);
      hPhi->Fill(phi);
      hTheta->Fill(theta);
      hImpulse->Fill(impulse);

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
        auto y = gRandom->Uniform(0, 1);
        Particle dau1{};
        Particle dau2{};
        if (0. <= y && y < 0.5) {
          dau1.setIndex("Pi+");
          dau2.setIndex("K-");
        } else {
          dau1.setIndex("Pi-");
          dau2.setIndex("K+");
        }
        if (eventParticles[i].Decay2Body(dau1, dau2) != 0) {
          break;
        } else {
          eventParticles[overflow].setIndex(dau1.getIndex());
          ++overflow;
          eventParticles[overflow].setIndex(dau2.getIndex());
          ++overflow;
          hInvMass5->Fill(dau1.particleInvMass(dau2));
        }
      }

      eventParticles[i].setImpulse(polarToCartesian(impulse, theta, phi));
      hTransvImpulse->Fill(
          sqrt(std::pow(impulse, 2) + std::pow(eventParticles[i].getImpulse().px_, 2)));

      hParticles->Fill(eventParticles[i].getIndex());

      hEnergy->Fill(eventParticles[i].particleEnergy());
    }

    for (int i{0}; i < overflow - 1; ++i) {
      for (int j{i + 1}; j < overflow; ++i) {
        hInvMass->Fill(eventParticles[i].particleInvMass(eventParticles[j]));
      }
    }

    for (int i{0}; i < overflow - 1; ++i) {
      if (eventParticles[i].getIndex() != 6) {
        for (int j{i + 1}; j < overflow; ++i) {
          if (eventParticles[j].getIndex() != 6) {
            if ((eventParticles[i].getIndex() + eventParticles[j].getIndex()) % 2 == 0)
              hInvMass1->Fill(eventParticles[i].particleInvMass(eventParticles[j]));
          } else {
            assert((eventParticles[i].getIndex() + eventParticles[j].getIndex()) % 2 == 1);
            hInvMass2->Fill(eventParticles[i].particleInvMass(eventParticles[j]));
          }
        }
      }
    }

    for (int i{0}; i < overflow - 1; ++i) {
      for (int j{i + 1}; j < overflow; ++j) {
        if ((eventParticles[i].getIndex() + eventParticles[j].getIndex()) == 3) {
          hInvMass3->Fill(eventParticles[i].particleInvMass(eventParticles[j]));
        } else if (eventParticles[i].getIndex() != eventParticles[j].getIndex()
                   && (((eventParticles[i].getIndex() + eventParticles[j].getIndex() == 2))
                       || ((eventParticles[i].getIndex() + eventParticles[j].getIndex()
                            == 4)))) {
          hInvMass4->Fill(eventParticles[i].particleInvMass(eventParticles[j]));
        }
      }
    }  
}
// TCanvas *myCanvas1 = new TCanvas("myCanvas1", "Dati: tipo particelle, angoli, impulso ed energia", 0, 10, 800, 600);
//     myCanvas1->Divide(2, 3); // dividi in due colonne con tre righe
//     for (int i{0}; i < 6; ++i)
//     {
//         myCanvas1->cd(i + 1);
//         hParticles->Draw();
//         hPhi->Draw()
//         hTheta->Draw()

//     }
//     // seconda canva per 6 histo
//     TCanvas *myCanvas2 = new TCanvas("myCanvas2", "Dati: massa invariante", 1000, 1000, 800, 600);
//     cHisto2->Divide(2, 3); // dividi in due colonne con tre righe
//     for (int j{6}; j < 12; ++j)
//     {
//         cHisto2->cd(j - 5);
//         h[j]->Draw("H");
//         h[j]->Draw("E,SAME");
//     }

hInvMass1->Draw();

TFile *file = new TFile("particle_simulation.root","RECREATE");
file->Write();
file->Close();
}