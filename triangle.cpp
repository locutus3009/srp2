#include <cmath>
#include "triangle.h"
#include "beam.h"

using namespace std;

// Triangle::Triangle(double x1_, double y1_, double z1_, 
// double x2_, double y2_, double z2_,
// double x3_, double y3_, double z3_) :  x1(x1_), y1(y1_), z1(z1_), x2(x2_), y2(y2_), z2(z2_), x3(x3_), y3(y3_), z3(z3_) {}
// {
//   SetTriangle(x1_, y1_, z1_, x2_, y2_, z2_, x3_, y3_, z3_);
//   a_Triangle(x1, y1, z1, x2, y2, z2);
//   b_Triangle(x1, y1, z1, x3, y3, z3);
//   squareTriangle = Get_squareTriangle(a1, a2, a3, b1, b2, b3);
//   Get_normalTriangle(a1, a2, a3, b1, b2, b3);
// }

void Triangle::Initialize(void)
{
  SetTriangle(x1, y1, z1, x2, y2, z2, x3, y3, z3);
  a_Triangle(x1, y1, z1, x2, y2, z2);
  b_Triangle(x1, y1, z1, x3, y3, z3);
  squareTriangle = Get_squareTriangle(a1, a2, a3, b1, b2, b3);
  Get_normalTriangle(a1, a2, a3, b1, b2, b3);
}

/*
*   строим треугольник
*/
void Triangle::SetTriangle(double x1_, double y1_, double z1_, 
double x2_, double y2_, double z2_,
double x3_, double y3_, double z3_)
{
  x1 = x1_; y1 = y1_; z1 = z1_;
  x2 = x2_; y2 = y2_; z2 = z2_;
  x3 = x3_; y3 = y3_; z3 = z3_;
}
/*
*   первый вектор
*/
void Triangle::a_Triangle(double x1_, double y1_, double z1_, 
double x2_, double y2_, double z2_)
{
  a1 = x2_ - x1_;
  a2 = y2_ - y1_;
  a3 = z2_ - z1_;
}
/*
*   второй вектор
*/
void Triangle::b_Triangle(double x1_, double y1_, double z1_, 
double x3_, double y3_, double z3_)
{
  b1 = x3_ - x1_;
  b2 = y3_ - y1_;
  b3 = z3_ - z1_;
}
/*
*   находим площадь треугольника
*/
double Triangle::Get_squareTriangle(double a1_, double a2_, double a3_, 
double b1_, double b2_, double b3_)
{
  return 0.5 * sqrt(pow((a2_*b3_ - a3_*b2_), 2) + pow((a3_*b1_ - a1_*b3_), 2) + pow((a1_*b2_ - a2_*b1_), 2));
}
/*
*   находим координаты нормали
*/
void Triangle::Get_normalTriangle(double a1_, double a2_, double a3_, 
double b1_, double b2_, double b3_)
{
  n1 = (a2_*b3_ - a3_*b2_) / sqrt(pow((a2_*b3_ - a3_*b2_), 2) + pow((a3_*b1_ - a1_*b3_), 2) + pow((a1_*b2_ - a2_*b1_), 2));
  n2 = (a3_*b1_ - a1_*b3_) / sqrt(pow((a2_*b3_ - a3_*b2_), 2) + pow((a3_*b1_ - a1_*b3_), 2) + pow((a1_*b2_ - a2_*b1_), 2));
  n3 = (a1_*b2_ - a2_*b1_) / sqrt(pow((a2_*b3_ - a3_*b2_), 2) + pow((a3_*b1_ - a1_*b3_), 2) + pow((a1_*b2_ - a2_*b1_), 2)); 

  // double norm = 1/sqrt(pow(n1,2)+pow(n2,2)+pow(n3,2));
  // n1 = n1*norm;
  // n2 = n2*norm;
  // n3 = n3*norm;
}
  
int Triangle::get_intersection(Beam beam, double *xp,double *yp,double *zp,double *t)
{
  //Константы в уравнении
  double A = y1*z2 - y1*z3 - y2*z1 + y2*z3 + y3*z1 - y3*z2;
  double B = -x1*z2 + x1*z3 + x2*z1 - x2*z3 - x3*z1 + x3*z2;
  double C = x1*y2 - x1*y3 - x2*y1 + x2*y3 + x3*y1 - x3*y2;
  double D = -x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1;
  //Расстояние
  double pogreshnost = MAX_ERROR;
  if(abs(A*beam.rx + B*beam.ry + C*beam.rz)<pogreshnost) return 0;
  *t = (-D - (A*beam.xq + B*beam.yq + C*beam.zq))/(A*beam.rx + B*beam.ry + C*beam.rz);
  if (abs(*t) < pogreshnost) return 0;
  if(*t < 0) return 0;
  //Координаты точки пересечения с плоскостью
  *xp = beam.xq + *t*beam.rx;
  *yp = beam.yq + *t*beam.ry;
  *zp = beam.zq + *t*beam.rz;
  
  //Проверка попадания в треугольник
  double xpa = x1 - *xp;
  double ypa = y1 - *yp;
  double zpa = z1 - *zp;
    
    
  double xpb = x2 - *xp;
  double ypb = y2 - *yp;
  double zpb = z2 - *zp;
  
  double xpc = x3 - *xp;
  double ypc = y3 - *yp;
  double zpc = z3 - *zp;
    
  double spapb = 0.5 * sqrt(pow((ypa*zpb - zpa*ypb), 2) + pow((zpa*xpb - xpa*zpb), 2) + pow((xpa*ypb - ypa*xpb), 2));
  double spapc = 0.5 * sqrt(pow((ypa*zpc - zpa*ypc), 2) + pow((zpa*xpc - xpa*zpc), 2) + pow((xpa*ypc - ypa*xpc), 2));
  double spbpc = 0.5 * sqrt(pow((ypb*zpc - zpb*ypc), 2) + pow((zpb*xpc - xpb*zpc), 2) + pow((xpb*ypc - ypb*xpc), 2));
  
  
  if (((spapb + spapc + spbpc) >= (squareTriangle - pogreshnost)) && ((spapb + spapc + spbpc) <= (squareTriangle + pogreshnost)))
  {
    return 1;
  }
  else return 0;
}
