#include "TRandom.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"

void loader(){
    gROOT->LoadMacro("particle_type.cpp+");
    gROOT->LoadMacro("resonance_type.cpp+");
    gROOT->LoadMacro("particle.cpp+");
}