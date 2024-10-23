#ifndef PARTICLE_TYPE_HPP
#define PARTICLE_TYPE_HPP
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <array>

class ParticleType
{
  const char *name_;
  const double mass_;
  const int charge_;

public:
  explicit ParticleType(char, double, int);
  const char *getName() const;
  const double getMass() const;
  const int getCharge() const;
  void virtual print() const;
};

#endif
