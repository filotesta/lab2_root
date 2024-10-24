#include "particle.hpp"

int Particle::nParticleType_ = 0;
std::array<ParticleType *, 10> Particle::ptrParticleType_{};

const int Particle::findParticle(const char *name)
{
  for (auto it : ptrParticleType_)
  {
    int index{0};
    if (*(it->getName()) == *name)
    {
      return index;
    }
    ++index;
  }
  return -1;
}

Particle::Particle(const char *name, Impulse impulse)
    : index_{findParticle(name)},
      impulse_{impulse}
{
}

Particle::Particle(const char *name)
    : Particle::Particle(name, Impulse{0., 0., 0.})
{
}

int Particle::getIndex()
{
  return index_;
}

void Particle::addParticleType(char *name, double mass, int charge, double width = 0)
{
  auto result{findParticle(name)};
  if (result == -1 && nParticleType_ < maxNumParticleType_)
  {
    if (width > 0)
    {
      ptrParticleType_[nParticleType_] = new ResonanceType{name, mass, charge, width};
    }
    else
    {
      ptrParticleType_[nParticleType_] = new ParticleType{name, mass, charge};
    }
    ++nParticleType_;
  }
  else
  {
    assert(*name == *(ptrParticleType_[result]->getName()));
  }
}

void Particle::setIndex(int index)
{
  if (index <= nParticleType_)
  {
    index_ = index;
  };
}

void Particle::setIndex(const char *name)
{
  auto index = Particle::findParticle(name);
  if (index <= nParticleType_)
  {
    index_ = index;
  };
}

const Impulse &Particle::getImpulse() const
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

double Particle::particleMass() const
{
  return ptrParticleType_[index_]->getMass();
};

double Particle::particleEnergy() const
{
  return sqrt(std::pow(particleMass(), 2) + std::pow(normImpulse(impulse_), 2));
}

double Particle::particleInvMass(const Particle &p) const
{
  return sqrt(std::pow((particleEnergy() + p.particleEnergy()), 2) - std::pow(normImpulse(sumVecImpulse(impulse_, p.getImpulse())), 2));
}