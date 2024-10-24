
#include "resonance_type.hpp"
#include <iostream>

ResonanceType::ResonanceType(const char* name, double mass, int charge, double width)
    : ParticleType(name, mass, charge)
    , width_{width}
{}

double ResonanceType::getWidth() const
{
  return width_;
}

void ResonanceType::print() const
{
  std::cout << "Resonance Type\n> Name: " << getName() << "\n> Mass: " << getMass()
            << "\n> Charge " << ResonanceType::getCharge()
            << "\n> Resonance width: " << width_ << '\n';
}
