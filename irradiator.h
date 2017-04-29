#ifndef irradiator_h_
#define irradiator_h_

#include <forward_list>
#include "beam.h"

class Irradiator
{
  public:
    std::forward_list<Beam> beams;

    Irradiator(void);
    void Init(Beam direction, double distance, double radius, int count);
    void Print(void);
};


#endif