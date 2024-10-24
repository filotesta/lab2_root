#ifndef RESONANCE_TYPE_HPP
#define RESONANCE_TYPE_HPP
#include "particle_type.hpp"

class ResonanceType : public ParticleType
{
  const double width_;

 public:
  explicit ResonanceType(const char*, double, int, double);
  double getWidth() const override;
  void print() const override;
};

#endif