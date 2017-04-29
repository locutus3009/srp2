#include "irradiator.h"
#include "beam.h"
#include <random>
#include <forward_list>
#include <cmath>
#include <chrono>

Irradiator::Irradiator(void)
{

}

void Irradiator::Init(Beam direction, double distance, double radius, int count)
{
  int i;
  beams.clear();
  Beam x(0, 0, 0, 0, 0, 1, 1,0);
  Beam y(0, 0, 0, 0, 0, 1, 1,0);

  if(direction.rz==0.0)
  {
    x.rx = 0.0;
    x.ry = 0.0;
    x.rz = 1.0;
  }
  else
  {
    x.rx = 0.0;
    x.ry = direction.rz/sqrt(pow(direction.ry,2)+pow(direction.rz,2));
    x.rz = -direction.ry/sqrt(pow(direction.ry,2)+pow(direction.rz,2));
  }

  y.rx = direction.ry*x.rz - direction.rz*x.ry;
  y.ry = direction.rz*x.rx - direction.rx*x.rz;
  y.rz = direction.rx*x.ry - direction.ry*x.rx;

  // x.Normalize();
  // y.Normalize();

  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 generator(seed);
  std::uniform_real_distribution<double> distribution(-1,1);

  double scale = 4*pow(radius,2)/count;

  auto it = beams.before_begin();
  for(i=0;i<count;++i)
  {
    double x_value = distribution(generator);
    double y_value = distribution(generator);
    double px = -distance*direction.rx + radius*x_value*x.rx + radius*y_value*y.rx;
    double py = -distance*direction.ry + radius*x_value*x.ry + radius*y_value*y.ry;
    double pz = -distance*direction.rz + radius*x_value*x.rz + radius*y_value*y.rz;
    it = beams.emplace_after ( it, px, py, pz, direction.rx, direction.ry, direction.rz, scale, 0);
  }
}

void Irradiator::Print(void)
{
  int i = 0;
  for ( auto it = beams.begin(); it != beams.end(); ++it )
  {
    printf("%d: %f %f %f %f %f %f %f %d\n",i,it->xq,it->yq,it->zq,it->rx,it->ry,it->rz,it->energy,it->counter);
    ++i;
  }
}