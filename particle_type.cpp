#include "particle_type.hpp"
#include <iostream>

ParticleType::ParticleType(const char *name, double mass, int charge)
    : name_{name}, mass_{mass}, charge_{charge}
{
}

const char *ParticleType::getName() const
{
  return name_;
}

double ParticleType::getMass() const
{
  return mass_;
}
int ParticleType::getCharge() const
{
  return charge_;
}

void ParticleType::print() const
{
  std::cout << "Particle Type\n> Name: " << name_ << "\n> Mass: " << mass_
            << "\n> Charge " << charge_ << "\n";
}
