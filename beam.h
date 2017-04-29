#ifndef beam_h_
#define beam_h_

#include <cmath>

class Beam
{
  public:

    double xq, yq, zq;
    double rx, ry, rz;
    double energy;
    int counter;
    /*
    *   Луч - координаты точки испускания и координаты направляющего вектора
    */
    Beam(double xq_, double yq_, double zq_, double rx_, double ry_, double rz_, double energy_, int counter_) :  xq(xq_), yq(yq_), zq(zq_), rx(rx_), ry(ry_), rz(rz_), energy(energy_), counter(counter_)
    {
          counter = counter_;
  //Нормализация
  double norm = 1/sqrt(pow(rx,2)+pow(ry,2)+pow(rz,2));
  rx = rx*norm;
  ry = ry*norm;
  rz = rz*norm;
    };
    void Normalize(void);
};

#endif