#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "resonance_type.hpp"
#include <array>
#include <cmath>
struct Impulse
{
  double px_;
  double py_;
  double pz_;
};

inline const double normImpulse(Impulse impulse)
{
  return sqrt(impulse.px_ * impulse.px_ + impulse.py_ * impulse.py_ + impulse.pz_ * impulse.pz_);
}

inline const Impulse& sumVecImpulse(const Impulse &p1, const Impulse &p2)
{
  return Impulse{p1.px_ + p2.px_, p1.py_ + p2.py_, p1.pz_ + p2.pz_};
}

class Particle
{
  static const int maxNumParticleType_{10};
  static std::array<ParticleType *, maxNumParticleType_> ptrParticleType_;
  static int nParticleType_;
  Impulse impulse_;
  int index_;
  static int findParticle(char name);

public:
  Particle(char, Impulse);
  const int getIndex();
  void setIndex(int);
  void setIndex(char); 
  const Impulse& getImpulse() const;
  void setImpulse(Impulse);
  static void addParticleType(char, double, int, double);
  static void printParticleTypes();
  void printParticleData();
  const double particleMass() const;
  const double particleEnergy() const;
  const double particleInvMass(const Particle &) const;
};

#endif