#include "particle_type.hpp"
#include <iostream>

explicit ParticleType::ParticleType(char name, double mass, int charge)
    : name_{&name}
    , mass_{mass}
    , charge_{charge}
{}

const char* ParticleType::getName() const
{
  return name_;
}

const double ParticleType::getMass() const
{
  return mass_;
}

const int ParticleType::getCharge() const
{
  return charge_;
}

void ParticleType::print() const
{
  std::cout << "Particle Type\n> Name: " << *name_ << "\n> Mass: " << mass_
            << "\n> Charge " << charge_ << "\n";
}
