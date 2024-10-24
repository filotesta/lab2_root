#include "particle.hpp"

int Particle::nParticleType_ = 0;
std::array<ParticleType*, Particle::maxNumParticleType_> Particle::ptrParticleType_{};

int Particle::findParticle(const char* name)
{
  for (int it{}; it < nParticleType_;) {
    if (ptrParticleType_[it]->getName() == name) {
      return it;
    }
    it++;
  }
  return -1;
}

Particle::Particle(const char* name, Impulse impulse)
    : index_{findParticle(name)}
    , impulse_{impulse}
{}

Particle::Particle(const char* name)
    : Particle::Particle(name, Impulse{0., 0., 0.})
{}

int Particle::getIndex()
{
  return index_;
}

void Particle::addParticleType(char* name, double mass, int charge, double width = 0)
{
  auto result{findParticle(name)};
  if (result == -1 && nParticleType_ < maxNumParticleType_) {
    if (width > 0) {
      ptrParticleType_[nParticleType_] = new ResonanceType{name, mass, charge, width};
    } else {
      ptrParticleType_[nParticleType_] = new ParticleType{name, mass, charge};
    }
    ++nParticleType_;
  } else {
    assert(*name == *(ptrParticleType_[result]->getName()));
  }
}

void Particle::setIndex(int index)
{
  if (index <= nParticleType_) {
    index_ = index;
  };
}

void Particle::setIndex(const char* name)
{
  auto index = Particle::findParticle(name);
  if (index <= nParticleType_) {
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

double Particle::particleMass() const
{
  return ptrParticleType_[index_]->getMass();
};

double Particle::particleEnergy() const
{
  return sqrt(std::pow(particleMass(), 2) + std::pow(normImpulse(impulse_), 2));
}

double Particle::particleInvMass(const Particle& p) const
{
  return sqrt(std::pow((particleEnergy() + p.particleEnergy()), 2)
              - std::pow(normImpulse(sumVecImpulse(impulse_, p.getImpulse())), 2));
} 

int Particle::Decay2Body(Particle& dau1, Particle& dau2) const
{
  if (particleMass() == 0.0) {
    printf("Decayment cannot be preformed if mass is zero\n");
    return 1;
  }
  double massMot  = particleMass();
  double massDau1 = dau1.particleMass();
  double massDau2 = dau2.particleMass();

  if (index_ > -1) { // add width effect

    // gaussian random numbers

    float x1, x2, w, y1;

    double invnum = 1. / RAND_MAX;
    do {
      x1 = 2.0 * rand() * invnum - 1.0;
      x2 = 2.0 * rand() * invnum - 1.0;
      w  = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w  = sqrt((-2.0 * log(w)) / w);
    y1 = x1 * w;

    massMot += ptrParticleType_[index_]->getWidth()
             * y1; // metodo getWidth anche in Partycle_type?
  }

  if (massMot < massDau1 + massDau2) {
    printf("Decayment cannot be preformed because mass is too low in this channel\n");
    return 2;
  }

  double pout =
      sqrt((massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2))
           * (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2)))
      / massMot * 0.5;

  double norm = 2 * M_PI / RAND_MAX;

  double phi   = rand() * norm;
  double theta = rand() * norm * 0.5 - M_PI / 2.;
  Impulse dau1Imp{pout * sin(theta) * cos(phi), pout * sin(theta) * sin(phi),
                  pout * cos(theta)};
  Impulse dau2Imp{-pout * sin(theta) * cos(phi), -pout * sin(theta) * sin(phi),
                  -pout * cos(theta)};
  dau1.setImpulse(dau1Imp);
  dau2.setImpulse(dau2Imp);

  double energy = sqrt(impulse_.px_ * impulse_.px_ + impulse_.py_ * impulse_.py_
                       + impulse_.pz_ * impulse_.pz_ + massMot * massMot);

  double bx = impulse_.px_ / energy;
  double by = impulse_.py_ / energy;
  double bz = impulse_.pz_ / energy;

  dau1.Boost(bx, by, bz);
  dau2.Boost(bx, by, bz);

  return 0;
}

void Particle::Boost(double bx, double by, double bz)
{
  double energy = particleEnergy();

  // Boost this Lorentz vector
  double b2     = bx * bx + by * by + bz * bz;
  double gamma  = 1.0 / sqrt(1.0 - b2);
  double bp     = bx * impulse_.px_ + by * impulse_.py_ + bz * impulse_.pz_;
  double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

  impulse_.px_ += gamma2 * bp * bx + gamma * bx * energy;
  impulse_.py_ += gamma2 * bp * by + gamma * by * energy;
  impulse_.pz_ += gamma2 * bp * bz + gamma * bz * energy;
}