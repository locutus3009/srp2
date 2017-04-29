#ifndef shape_h_
#define shape_h_

#include <forward_list>
#include "triangle.h"
#include "beam.h"

class Shape
{
  public:
    int count;
    std::forward_list<Triangle> triangles;

    Shape(void);
    void Load(char* filename);
    void Print(void);
    Triangle * get_intersection(Beam ray,double * t,double * Px,double * Py,double * Pz);
};

#endif