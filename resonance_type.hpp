#ifndef RESONANCE_TYPE_HPP
#define RESONANCE_TYPE_HPP
#include "particle_type.hpp"

class ResonanceType : public ParticleType
{
  const double width_;

 public:
  explicit ResonanceType(char name, double mass, int charge, double width);
  const double getWidth() const;
  void print() const override;
};


#endif