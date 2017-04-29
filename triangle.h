#ifndef triangle_h_
#define triangle_h_

#include "beam.h"

#define MAX_ERROR 0.000001




class Triangle
{
  public:

    //Вершины треугольника
    double x1, y1, z1;
    double x2, y2, z2;
    double x3, y3, z3;
    //вектора
    double a1, a2, a3;
    double b1, b2, b3;
    //нормаль
    double n1, n2, n3;
    //Оптические параметры
    double rho,s,B;

    double squareTriangle;
    /*
    *   получаем координаты точек, строим треугольник, находим координаты векторов-сторон, находим площадь, находим нормаль,
    *   выводим площадь на монитор
    */
    Triangle(double x1_, double y1_, double z1_, 
        double x2_, double y2_, double z2_,
        double x3_, double y3_, double z3_) :  x1(x1_), y1(y1_), z1(z1_), x2(x2_), y2(y2_), z2(z2_), x3(x3_), y3(y3_), z3(z3_) {};
    void Initialize(void);
    void SetTriangle(double x1_, double y1_, double z1_, 
        double x2_, double y2_, double z2_,
        double x3_, double y3_, double z3_);
    void a_Triangle(double x1_, double y1_, double z1_, 
        double x2_, double y2_, double z2_);
    void b_Triangle(double x1_, double y1_, double z1_, 
        double x3_, double y3_, double z3_);
    double Get_squareTriangle(double a1_, double a2_, double a3_, 
        double b1_, double b2_, double b3_);
    void Get_normalTriangle(double a1_, double a2_, double a3_, 
        double b1_, double b2_, double b3_);
    int get_intersection(Beam beam, double *xp,double *yp,double *zp,double *t);
};

#endif