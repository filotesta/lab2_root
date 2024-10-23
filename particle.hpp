#ifndef PARTICLE_HPP
#define PARTICLE_HPP
#include "resonance_type.hpp"

struct Impulse
{
  double px_;
  double py_;
  double pz_;
};

inline double normImpulse(Impulse impulse)
{
  return sqrt(impulse.px_ * impulse.px_ + impulse.py_ * impulse.py_ + impulse.pz_ * impulse.pz_);
}

inline Impulse sumVecImpulse(const Impulse &p1, const Impulse &p2)
{
  return Impulse{p1.px_ + p2.px_, p1.py_ + p2.py_, p1.pz_ + p2.pz_};
}

class Particle
{
  static const int maxNumParticleType_{10};
  static std::array<ParticleType *, maxNumParticleType_> ptrParticleType_;
  static int nParticleType_;
  int index_;
  Impulse impulse_;
  static int findParticle(char name);

public:
  Particle(char);
  Particle(char, Impulse);
  int getIndex();
  void setIndex(int);
  void setIndex(char);
  const Impulse &getImpulse() const;
  void setImpulse(Impulse);
  static void addParticleType(char, double, int, double);
  static void printParticleTypes();
  void printParticleData();
  double particleMass() const;
  double particleEnergy() const;
  double particleInvMass(const Particle &) const;
};

#endif