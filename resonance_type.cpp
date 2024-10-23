
#include "resonance_type.hpp"
#include <iostream>

ResonanceType::ResonanceType(char name, double mass, int charge, double width)
    : ParticleType(name, mass, charge)
    , width_{width}
{}

double ResonanceType::getWidth() const
{
  return width_;
}

void ResonanceType::print() const
{
  ParticleType::print();
  std::cout << "> Resonance width: " << width_ << "\n";
}

