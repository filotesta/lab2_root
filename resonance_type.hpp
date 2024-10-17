#ifndef RESONANCE_TYPE.HPP
#define RESONANCE_TYPE .HPP

#include "particle_type.hpp"
class ResonanceType : public ParticleType
{
  const double width_;

 public:
  ResonanceType(char name, double mass, int charge, double width);
  const double getWidth() const;
  void print() const override;
};


#endif