#ifndef PARTICLE_TYPE_HPP
#define PARTICLE_TYPE_HPP
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>

class ParticleType
{
  const char* name_;
  const double mass_;
  const int charge_;

 public:
  explicit ParticleType(const char*, double, int);
  const char* getName() const;
  double getMass() const;
  int getCharge() const;
  void virtual print() const;
  double virtual getWidth() const;
};

#endif
