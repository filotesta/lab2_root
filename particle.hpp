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
  return sqrt(impulse.px_ * impulse.px_ + impulse.py_ * impulse.py_
              + impulse.pz_ * impulse.pz_);
}

inline Impulse sumVecImpulse(const Impulse& p1, const Impulse& p2)
{
  return Impulse{p1.px_ + p2.px_, p1.py_ + p2.py_, p1.pz_ + p2.pz_};
}

class Particle
{
 public:
  explicit Particle();
  explicit Particle(const char*);
  explicit Particle(const char*, Impulse);
  static const int maxNumParticleType_{10};
  int getIndex();
  void setIndex(int);
  void setIndex(const char*);
  const Impulse& getImpulse() const;
  void setImpulse(Impulse);
  static void addParticleType(const char*, double, int, double);
  static void printParticleTypes();
  void printParticleData();
  double particleMass() const;
  double particleEnergy() const;
  double particleInvMass(const Particle&) const;
  int getCharge() const;

  int Decay2Body(Particle&, Particle&) const;
  void Boost(double bx, double by, double bz);

 private:
  static std::array<ParticleType*, maxNumParticleType_> ptrParticleType_;
  static int nParticleType_;
  int index_;
  Impulse impulse_;
  static int findParticle(const char* name);
};

#endif