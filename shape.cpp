#include "shape.h"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <forward_list>
#include "beam.h"
#include "triangle.h"

using namespace std;

Shape::Shape(void)
{
  count = 0;
}

void Shape::Load(char* filename)
{
  if(filename==NULL) return;
  FILE *in = fopen( (const char*) filename, "r" );
  triangles.clear();
  auto it = triangles.before_begin();
  while(!feof(in))
  {
    float ax,ay,az,bx,by,bz,cx,cy,cz,rho,s,B;
    fscanf(in,"%f %f %f %f %f %f %f %f %f %f %f %f",&ax,&ay,&az,&bx,&by,&bz,&cx,&cy,&cz,&rho,&s,&B);
    if(feof(in)) return;
    // printf("%d: %f %f %f %f %f %f %f %f %f\n",count,ax, ay, az, bx, by, bz, cx, cy, cz);
    it = triangles.emplace_after ( it, ax, ay, az, bx, by, bz, cx, cy, cz);
    it->Initialize();
    it->rho = rho;
    it->s = s;
    it->B = B;
    ++count;
  }
  fclose(in);
}

void Shape::Print(void)
{
  int i = 0;
  for ( auto it = triangles.begin(); it != triangles.end(); ++it )
  {
    printf("%d: %f %f %f %f %f %f %f %f %f %f %f %f %f\n",i,it->squareTriangle,it->x1,it->y1,it->z1,it->x2,it->y2,it->z2,it->x3,it->y3,it->z3,it->rho,it->s,it->B);
    ++i;
  }

}

Triangle * Shape::get_intersection(Beam ray,double * t,double * Px,double * Py,double * Pz)
{
  forward_list< pair<Triangle, pair<double, pair<double, pair<double, double > > > > > selected_triangles;
  auto selected_triangles_iterator = selected_triangles.before_begin();


  for ( auto it = triangles.begin(); it != triangles.end(); ++it )
  {
    double t,xp,yp,zp;
    if(it->get_intersection(ray,&xp,&yp,&zp,&t) == 1)
    {
      selected_triangles_iterator = selected_triangles.insert_after ( selected_triangles_iterator, 1, {*it,{t,{xp,{yp,zp}}}});
    }
  }
  if (selected_triangles.empty()) return NULL;

  int i = 0;
  double t_new = 1000000.0;
  double t_old = t_new;
  Triangle * tr_new;
  for ( auto it = selected_triangles.begin(); it != selected_triangles.end(); ++it )
  {
    // printf("%d: %f %f %f %f %f %f %f %f %f %f %f %f %f\n",i,it->first.squareTriangle,it->first.x1,it->first.y1,it->first.z1,it->first.x2,it->first.y2,it->first.z2,it->first.x3,it->first.y3,it->first.z3,it->first.rho,it->first.s,it->first.B);
    // printf("t: %f, Px: %f, Py: %f, Pz: %f\n",it->second.first,it->second.second.first,it->second.second.second.first,it->second.second.second.second);
    ++i;
    t_new = it->second.first;
    if(t_new<t_old)
    {
      t_old = t_new;
      tr_new = &(it->first);
      *t = t_new;
      *Px = it->second.second.first;
      *Py = it->second.second.second.first;
      *Pz = it->second.second.second.second;
    }
  }
  return tr_new;
}