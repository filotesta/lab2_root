#include "particle.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

int Particle::nParticleType_ = 0;

int Particle::findParticle(char name)
{
  for (auto it : ptrParticleType_)
  {
    int index{0};
    if (*(it->getName()) == name)
    {
      return index;
    }
    ++index;
  }
  return -1;
}
Particle::Particle(char name, Impulse impulse = {0, 0, 0})
    : index_{findParticle(name)}, impulse_{impulse}
{
}

const int Particle::getIndex()
{
  return index_;
}

void Particle::addParticleType(char name, double mass, int charge, double width = 0)
{
  auto result{findParticle(name)};
  if (result == -1 && nParticleType_ < maxNumParticleType_)
  {
    if (width > 0)
    {
      ptrParticleType_[nParticleType_] = &ResonanceType(name, mass, charge, width);
    }
    else
    {
      ptrParticleType_[nParticleType_] = &ParticleType(name, mass, charge);
    }
    ++nParticleType_;
  }
  else
  {
    assert(name == *(ptrParticleType_[result]->getName()));
  }
}

void Particle::setIndex(int index)
{
  if (index <= nParticleType_)
  {
    index_ = index;
  };
}

void Particle::setIndex(char name)
{
  auto index = Particle::findParticle(name);
  if (index <= nParticleType_)
  {
    index_ = index;
  };
}

const Impulse& Particle::getImpulse() const
{
  return impulse_;
}

void Particle::setImpulse(Impulse impulse)
{
  impulse_ = impulse;
}

void Particle::printParticleTypes()
{
  for (auto it : ptrParticleType_)
    it->print();
}

void Particle::printParticleData()
{
  std::cout << "particle's name: " << ptrParticleType_[index_]->getName()
            << "\nindex: " << index_ << "\nimpulse: ( " << this->impulse_.px_ << ", "
            << this->impulse_.py_ << ", " << this->impulse_.pz_ << " )\n";
}

const double Particle::particleMass() const
{
  return ptrParticleType_[index_]->getMass();
};

const double Particle::particleEnergy() const
{
  return std::sqrt(std::pow(particleMass(), 2) + std::pow(normImpulse(impulse_), 2));
}

const double Particle::particleInvMass(const Particle &p) const
{
  return std::sqrt(std::pow((particleEnergy() + p.particleEnergy()), 2) - std::pow(normImpulse(sumVecImpulse(impulse_, p.getImpulse())), 2));
}