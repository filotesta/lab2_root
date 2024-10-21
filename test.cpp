#include "particle.hpp"
#include <iostream>

int main()
{
    ParticleType Kplus{'K+', 0.500, -2};
    ParticleType Kminus{'K-', 0.400, -3};
    ParticleType Pminus{'P-', 0.400, -3};
    ParticleType Pplus{'P-', 0.400, -3};
    ResonanceType Pplus2{'P-', 0.400, -3, 5.45};
    Kplus.print();
    Pplus2.print();
    Impulse impulse{0, 0.34, 0.56};
    Particle kkplus{'K+',impulse};
    Particle kkminus{'K-',impulse};
    std::cout<<kkminus.getIndex()<<'\n';    
    kkminus.setIndex(kkplus.getIndex());
    std::cout<<kkminus.getIndex()<<'\n';
}