#include <cmath>
#include "beam.h"

using namespace std;

// Beam::Beam(double xq_, double yq_, double zq_, double rx_, double ry_, double rz_)
// {
//   xq = xq_;
//   yq = yq_;
//   zq = zq_;
//   counter = 0;
//   //Нормализация
//   double norm = 1/sqrt(pow(rx_,2)+pow(ry_,2)+pow(rz_,2));
//   rx = rx_*norm;
//   ry = ry_*norm;
//   rz = rz_*norm;
// }

void Beam::Normalize(void)
{
  double norm = 1/sqrt(pow(rx,2)+pow(ry,2)+pow(rz,2));
  rx = rx*norm;
  ry = ry*norm;
  rz = rz*norm;
}