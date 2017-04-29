#include <cmath>
#include <cstdio>
#include <iostream>
#include <forward_list>
#include <random>
#include <chrono>
#include <cstring>
#include "beam.h"
#include "triangle.h"
#include "shape.h"
#include "irradiator.h"



#define TEST 0


using namespace std;


void test_1(void)
{
  double xp;
  double yp;
  double zp;
  double t;

  Triangle Tri_rod(0, 0, 0, 3, 0, 0, 0, 3, 0);//test triangle
  Tri_rod.Initialize();
  Beam Beam_rod(1, 1, 1, 0, 0, -1, 1, 0);
  //checking
  if (Tri_rod.get_intersection(Beam_rod,&xp,&yp,&zp,&t) == 1)
  {
    cout << "INSIDE\n";
    printf("(%f,%f,%f), t=%f\n",xp,yp,zp,t);
  }
  else
  {
    cout << "OUTSIDE\n";
  }
}

void test_2(void)
{
  double xp;
  double yp;
  double zp;
  double t;

  Triangle Tri_rod(0, 0, 0, 3, 0, 0, 0, 3, 0);//test triangle
  Tri_rod.Initialize();
  Beam Beam_rod(3, 1, 1, 0, 0, -1, 1, 0);
  //checking
  if (Tri_rod.get_intersection(Beam_rod,&xp,&yp,&zp,&t) == 1)
  {
    cout << "INSIDE\n";
    printf("(%f,%f,%f), t=%f\n",xp,yp,zp,t);
  }
  else
  {
    cout << "OUTSIDE\n";
  }
}

void test_3(void)
{
  double xp;
  double yp;
  double zp;
  double t;

  Triangle Tri_rod(0, 0, 0, 3, 0, 0, 0, 3, 0);//test triangle
  Tri_rod.Initialize();
  Beam Beam_rod(1, 0, 6, 0, 0, -1 ,1, 0);
  //checking
  if (Tri_rod.get_intersection(Beam_rod,&xp,&yp,&zp,&t) == 1)
  {
    cout << "INSIDE\n";
    printf("(%f,%f,%f), t=%f\n",xp,yp,zp,t);
  }
  else
  {
    cout << "OUTSIDE\n";
  }
}

void test_4(void)
{
  double xp;
  double yp;
  double zp;
  double t;

  Triangle Tri_rod(0, 0, 0, 3, 0, 0, 0, 3, 0);//test triangle
  Tri_rod.Initialize();
  Beam Beam_rod(2, 1, 1, -1, 0, -1, 1, 0);
  //checking
  if (Tri_rod.get_intersection(Beam_rod,&xp,&yp,&zp,&t) == 1)
  {
    cout << "INSIDE\n";
    printf("(%f,%f,%f), t=%f\n",xp,yp,zp,t);
  }
  else
  {
    cout << "OUTSIDE\n";
  }
}

void test_5(void)
{
  double xp;
  double yp;
  double zp;
  double t;

  Triangle Tri_rod(0, 0, 0, 3, 0, 1, 0, 3, 0);//test triangle
  Tri_rod.Initialize();
  Beam Beam_rod(2, 1, 1, -1, 0, -1, 1, 0);
  //checking
  if (Tri_rod.get_intersection(Beam_rod,&xp,&yp,&zp,&t) == 1)
  {
    cout << "INSIDE\n";
    printf("(%f,%f,%f), t=%f\n",xp,yp,zp,t);
  }
  else
  {
    cout << "OUTSIDE\n";
  }
}

void test_6(void)
{
  double xp;
  double yp;
  double zp;
  double t;

  Triangle tr1(0, 0, 0, 3, 0, 0, 0, 3, 0);//test triangle
  tr1.Initialize();
  Triangle tr2(0, 0, 1, 3, 0, 1, 0, 3, 1);//test triangle
  tr2.Initialize();
  Beam Beam_rod(1, 1, 3, 0, 0, -1, 1, 0);
  //checking
  if (tr1.get_intersection(Beam_rod,&xp,&yp,&zp,&t) == 1)
  {
    cout << "INSIDE\n";
    printf("(%f,%f,%f), t=%f\n",xp,yp,zp,t);
  }
  else
  {
    cout << "OUTSIDE\n";
  }
  if (tr2.get_intersection(Beam_rod,&xp,&yp,&zp,&t) == 1)
  {
    cout << "INSIDE\n";
    printf("(%f,%f,%f), t=%f\n",xp,yp,zp,t);
  }
  else
  {
    cout << "OUTSIDE\n";
  }
}

void test_7(void)
{
  double xp;
  double yp;
  double zp;
  double t;
  int i = 0;

  Triangle Tri_rod(0, 0, 0, 3, 0, 1, 0, 3, 0);//test triangle
  Tri_rod.Initialize();
  Beam Beam_rod(2, 1, 1, -1, 0, -1, 1, 0);

  forward_list<Beam> beams = {Beam_rod};

  for ( auto it = beams.begin(); it != beams.end(); ++it )
  {
    printf("%d: ",i);
    ++i;
    if (Tri_rod.get_intersection(*it,&xp,&yp,&zp,&t) == 1)
    {
      cout << "INSIDE\n";
    }
    else
    {
      cout << "OUTSIDE\n";
    }
    if(i<=3)
    {
      beams.insert_after ( it, 3, Beam_rod );
    }
  }

}

void test_8(void)
{
  Shape * shape = new Shape();
  shape->Load("../examples/plane4.in");
  shape->Print();
}

void test_9(void)
{
  printf("==========TEST 09==========\n");
  Shape * shape = new Shape();
  shape->Load("../examples/plane1_new.in");
  shape->Print();
  double t,Px,Py,Pz;
  Triangle * tri;
  Beam beam(0.4, 0.5, 1, -0, 0, -1, 1, 0);
  tri = shape->get_intersection(beam,&t,&Px,&Py,&Pz);
  if(tri!=NULL)
  {
    printf("Intersection found\nt: %f Px: %f Py: %f Pz: %f\n",t,Px,Py,Pz);
  }
}

void test_10(void)
{
  printf("==========TEST 10==========\n");
  Irradiator irradiator;
  irradiator.Init({0,0,0,0,0,-1,1,0},10,1,10);
  irradiator.Print();
}

void test_11(void)
{
  printf("==========TEST 11==========\nAbsorbing plane\n");
  Irradiator irradiator;
  irradiator.Init({0,0,0,0,0,-1,1,0},10,1.0,100);
  // irradiator.Print();
  Shape * shape = new Shape();
  shape->Load("../examples/plane1_new.in");

  double Fx = 0.0;
  double Fy = 0.0;
  double Fz = 0.0;

  Triangle * tri;
  double t,Px,Py,Pz;

  for ( auto it = irradiator.beams.begin(); it != irradiator.beams.end(); ++it )
  {
    tri = shape->get_intersection(*it,&t,&Px,&Py,&Pz);
    if(tri!=NULL)
    {
      // printf("Intersection found: t: %f Px: %f Py: %f Pz: %f\n",t,Px,Py,Pz);
      Fx = Fx + (it->energy)*(it->rx);
      Fy = Fy + (it->energy)*(it->ry);
      Fz = Fz + (it->energy)*(it->rz);
    }
  }
  printf("Fx: %f, Fy: %f, Fz: %f\n",Fx,Fy,Fz);
}

void test_12(void)
{
  printf("==========TEST 12==========\nSpecular plane\n");
  Irradiator irradiator;
  irradiator.Init({0,0,0,0,0,-1,1,0},10,2.0,100);
  // irradiator.Print();
  Shape * shape = new Shape();
  shape->Load("../examples/plane1_new.in");

  double Fx = 0.0;
  double Fy = 0.0;
  double Fz = 0.0;

  Triangle * tri;
  double t,Px,Py,Pz;

  for ( auto it = irradiator.beams.begin(); it != irradiator.beams.end(); ++it )
  {
    tri = shape->get_intersection(*it,&t,&Px,&Py,&Pz);
    if(tri!=NULL)
    {
      printf("Intersection found: t: %f Px: %f Py: %f Pz: %f\n",t,Px,Py,Pz);
      Fx = Fx + (it->energy)*(it->rx);
      Fy = Fy + (it->energy)*(it->ry);
      Fz = Fz + (it->energy)*(it->rz);
      double nr1 = (tri->n1)*(it->rx) + (tri->n2)*(it->ry) + (tri->n3)*(it->rz);
      double sign = 1.0;
      if(nr1>0) sign = -1.0;
      double r2x = it->rx - sign*2*nr1*(tri->n1);
      double r2y = it->ry - sign*2*nr1*(tri->n2);
      double r2z = it->rz - sign*2*nr1*(tri->n3);
      if((it->counter)<10)
      {
        auto itnew = irradiator.beams.emplace_after ( it, Px, Py, Pz,r2x, r2y, r2z, (it->energy)*(tri->s)*(tri->rho), (it->counter)+1);
        Fx = Fx - (itnew->energy)*(itnew->rx);
        Fy = Fy - (itnew->energy)*(itnew->ry);
        Fz = Fz - (itnew->energy)*(itnew->rz);
      }
    }
  }
  printf("Fx: %f, Fy: %f, Fz: %f\n",Fx,Fy,Fz);
}

void test_13(void)
{
  printf("==========TEST 13==========\nDiffuse plane\n");
  Irradiator irradiator;
  irradiator.Init({0,0,0,0,0,-1,1,0},10,2.0,10000);
  // irradiator.Print();
  Shape * shape = new Shape();
  shape->Load("../examples/plane10_diffuse.in");

  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 generator(seed);
  std::uniform_real_distribution<double> distribution(-1,1);

  double Fx = 0.0;
  double Fy = 0.0;
  double Fz = 0.0;

  Triangle * tri;
  double t,Px,Py,Pz;

  for ( auto it = irradiator.beams.begin(); it != irradiator.beams.end(); ++it )
  {
    tri = shape->get_intersection(*it,&t,&Px,&Py,&Pz);
    if(tri!=NULL)
    {
      Fx = Fx + (it->energy)*(it->rx);
      Fy = Fy + (it->energy)*(it->ry);
      Fz = Fz + (it->energy)*(it->rz);
      double nr1 = (tri->n1)*(it->rx) + (tri->n2)*(it->ry) + (tri->n3)*(it->rz);
      double r2x = it->rx - 2*nr1*(tri->n1);
      double r2y = it->ry - 2*nr1*(tri->n2);
      double r2z = it->rz - 2*nr1*(tri->n3);

      double r3x = distribution(generator);
      double r3y = distribution(generator);
      double r3z_a = distribution(generator);
      double r3z = r3z_a/abs(r3z_a)*sqrt(pow(r3x,2) + pow(r3y,2));
      while((r3x*it->rx + r3y*it->ry + r3z*it->rz)>0)
      {
        r3x = distribution(generator);
        r3y = distribution(generator);
        r3z_a = distribution(generator);
        r3z = r3z_a/abs(r3z_a)*sqrt(pow(r3x,2) + pow(r3y,2));
      }

      if((it->counter)<10)
      {
        auto itnew = irradiator.beams.emplace_after ( it, Px, Py, Pz,r2x, r2y, r2z, (it->energy)*(tri->s)*(tri->rho), (it->counter)+1);
        Fx = Fx - (itnew->energy)*(itnew->rx);
        Fy = Fy - (itnew->energy)*(itnew->ry);
        Fz = Fz - (itnew->energy)*(itnew->rz);
      }
      if((it->counter)<10)
      {
        auto itnew = irradiator.beams.emplace_after ( it, Px, Py, Pz,r3x, r3y, r3z, (it->energy)*(1 - tri->s)*(tri->rho), (it->counter)+1);
        Fx = Fx - (itnew->energy)*(itnew->rx);
        Fy = Fy - (itnew->energy)*(itnew->ry);
        Fz = Fz - (itnew->energy)*(itnew->rz);
      }
    }
  }
  printf("Fx: %f, Fy: %f, Fz: %f\n",Fx,Fy,Fz);
}

void test_14(void)
{
  printf("==========TEST 14==========\nSpecular-diffuse plane\n");
  Irradiator irradiator;
  irradiator.Init({0,0,0,0,0,1,1,0},10,2.0,10000);
  // irradiator.Print();
  Shape * shape = new Shape();
  shape->Load("../examples/plane10_diffuse_specular.in");

  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 generator(seed);
  std::uniform_real_distribution<double> distribution(-1,1);

  double Fx = 0.0;
  double Fy = 0.0;
  double Fz = 0.0;

  Triangle * tri;
  double t,Px,Py,Pz;

  for ( auto it = irradiator.beams.begin(); it != irradiator.beams.end(); ++it )
  {
    tri = shape->get_intersection(*it,&t,&Px,&Py,&Pz);
    if(tri!=NULL)
    {
      Fx = Fx + (it->energy)*(it->rx);
      Fy = Fy + (it->energy)*(it->ry);
      Fz = Fz + (it->energy)*(it->rz);
      double nr1 = (tri->n1)*(it->rx) + (tri->n2)*(it->ry) + (tri->n3)*(it->rz);
      double r2x = it->rx - 2*nr1*(tri->n1);
      double r2y = it->ry - 2*nr1*(tri->n2);
      double r2z = it->rz - 2*nr1*(tri->n3);

      double r3x = distribution(generator);
      double r3y = distribution(generator);
      double r3z_a = distribution(generator);
      double r3z = r3z_a/abs(r3z_a)*sqrt(pow(r3x,2) + pow(r3y,2));
      while((r3x*it->rx + r3y*it->ry + r3z*it->rz)>0)
      {
        r3x = distribution(generator);
        r3y = distribution(generator);
        r3z_a = distribution(generator);
        r3z = r3z_a/abs(r3z_a)*sqrt(pow(r3x,2) + pow(r3y,2));
      }

      if((it->counter)<10)
      {
        auto itnew = irradiator.beams.emplace_after ( it, Px, Py, Pz,r2x, r2y, r2z, (it->energy)*(tri->s)*(tri->rho), (it->counter)+1);
        Fx = Fx - (itnew->energy)*(itnew->rx);
        Fy = Fy - (itnew->energy)*(itnew->ry);
        Fz = Fz - (itnew->energy)*(itnew->rz);
      }
      if((it->counter)<10)
      {
        auto itnew = irradiator.beams.emplace_after ( it, Px, Py, Pz,r3x, r3y, r3z, (it->energy)*(1 - tri->s)*(tri->rho), (it->counter)+1);
        Fx = Fx - (itnew->energy)*(itnew->rx);
        Fy = Fy - (itnew->energy)*(itnew->ry);
        Fz = Fz - (itnew->energy)*(itnew->rz);
      }
    }
  }
  printf("Fx: %f, Fy: %f, Fz: %f\n",Fx,Fy,Fz);
}




void trace(char * input,double dir_x, double dir_y, double dir_z,double distance,double radius,int count,char savetofile)
{
  Irradiator irradiator;
  irradiator.Init({0,0,0,dir_x,dir_y,dir_z,1,0},distance,radius,count);
  Shape * shape = new Shape();
  shape->Load(input);

  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 generator(seed);
  std::uniform_real_distribution<double> distribution(-1,1);

  double Fx = 0.0;
  double Fy = 0.0;
  double Fz = 0.0;

  Triangle * tri;
  double t,Px,Py,Pz;

  for ( auto it = irradiator.beams.begin(); it != irradiator.beams.end(); ++it )
  {
    tri = shape->get_intersection(*it,&t,&Px,&Py,&Pz);
    if(tri!=NULL)
    {
      Fx = Fx + (it->energy)*(it->rx);
      Fy = Fy + (it->energy)*(it->ry);
      Fz = Fz + (it->energy)*(it->rz);
      double nr1 = (tri->n1)*(it->rx) + (tri->n2)*(it->ry) + (tri->n3)*(it->rz);
      double r2x = it->rx - 2*nr1*(tri->n1);
      double r2y = it->ry - 2*nr1*(tri->n2);
      double r2z = it->rz - 2*nr1*(tri->n3);

      if(((it->counter)<10)&&((tri->s)>0.0))
      {
        auto itnew = irradiator.beams.emplace_after ( it, Px, Py, Pz,r2x, r2y, r2z, (it->energy)*(tri->s)*(tri->rho), (it->counter)+1);
        Fx = Fx - (itnew->energy)*(itnew->rx);
        Fy = Fy - (itnew->energy)*(itnew->ry);
        Fz = Fz - (itnew->energy)*(itnew->rz);
      }
      for(int j=0;j<1;++j)
      {
        double r3x = distribution(generator);
        double r3y = distribution(generator);
        double r3z_a = distribution(generator);
        while(r3z_a==0.0)
        {
          r3z_a = distribution(generator);
        }
        double r3z = r3z_a/abs(r3z_a)*sqrt(pow(r3x,2) + pow(r3y,2));
        while((r3x*it->rx + r3y*it->ry + r3z*it->rz)>0)
        {
          r3x = distribution(generator);
          r3y = distribution(generator);
          r3z_a = distribution(generator);
          while(r3z_a==0.0)
          {
            r3z_a = distribution(generator);
          }
          r3z = r3z_a/abs(r3z_a)*sqrt(pow(r3x,2) + pow(r3y,2));
        }
        if(((it->counter)<10)&&((tri->s)<1.0))
        {
          auto itnew = irradiator.beams.emplace_after ( it, Px, Py, Pz,r3x, r3y, r3z, (it->energy)*(1 - tri->s)*(tri->rho), (it->counter)+1);
          Fx = Fx - (itnew->energy)*(itnew->rx);
          Fy = Fy - (itnew->energy)*(itnew->ry);
          Fz = Fz - (itnew->energy)*(itnew->rz);
        }
      }
    }
  }
  if(savetofile==0)
  {
    printf("%f %f %f\n",Fx,Fy,Fz);
  }
  else
  {
    printf("%f %f %f\n",Fx,Fy,Fz);
    FILE *out = fopen( "output.csv", "a" );
    fprintf(out,"%f,%f,%f\n",Fx,Fy,Fz);
    fclose(out);
  }
}



void trace2(Shape * shape, double theta, double beta, double distance,double radius,int count,char savetofile,FILE *out)
{
  double dir_x = -sin(beta)*cos(theta);
  double dir_y = -sin(beta)*sin(theta);
  double dir_z = -cos(beta);

  Irradiator irradiator;
  irradiator.Init({0,0,0,dir_x,dir_y,dir_z,1,0},distance,radius,count);


  auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::mt19937 generator(seed);
  std::uniform_real_distribution<double> distribution(-1,1);

  double Fx = 0.0;
  double Fy = 0.0;
  double Fz = 0.0;

  Triangle * tri;
  double t,Px,Py,Pz;

  for ( auto it = irradiator.beams.begin(); it != irradiator.beams.end(); ++it )
  {
    tri = shape->get_intersection(*it,&t,&Px,&Py,&Pz);
    if(tri!=NULL)
    {
      Fx = Fx + (it->energy)*(it->rx);
      Fy = Fy + (it->energy)*(it->ry);
      Fz = Fz + (it->energy)*(it->rz);
      double nr1 = (tri->n1)*(it->rx) + (tri->n2)*(it->ry) + (tri->n3)*(it->rz);
      double r2x = it->rx - 2*nr1*(tri->n1);
      double r2y = it->ry - 2*nr1*(tri->n2);
      double r2z = it->rz - 2*nr1*(tri->n3);

      if(((it->counter)<10)&&((tri->s)>0.0))
      {
        auto itnew = irradiator.beams.emplace_after ( it, Px, Py, Pz,r2x, r2y, r2z, (it->energy)*(tri->s)*(tri->rho), (it->counter)+1);
        Fx = Fx - (itnew->energy)*(itnew->rx);
        Fy = Fy - (itnew->energy)*(itnew->ry);
        Fz = Fz - (itnew->energy)*(itnew->rz);
      }
      for(int j=0;j<1;++j)
      {
        double r3x = distribution(generator);
        double r3y = distribution(generator);
        double r3z_a = distribution(generator);
        while(r3z_a==0.0)
        {
          r3z_a = distribution(generator);
        }
        double r3z = r3z_a/abs(r3z_a)*sqrt(pow(r3x,2) + pow(r3y,2));
        while((r3x*it->rx + r3y*it->ry + r3z*it->rz)>0)
        {
          r3x = distribution(generator);
          r3y = distribution(generator);
          r3z_a = distribution(generator);
          while(r3z_a==0.0)
          {
            r3z_a = distribution(generator);
          }
          r3z = r3z_a/abs(r3z_a)*sqrt(pow(r3x,2) + pow(r3y,2));
        }
        if(((it->counter)<10)&&((tri->s)<1.0))
        {
          auto itnew = irradiator.beams.emplace_after ( it, Px, Py, Pz,r3x, r3y, r3z, (it->energy)*(1 - tri->s)*(tri->rho), (it->counter)+1);
          Fx = Fx - (itnew->energy)*(itnew->rx);
          Fy = Fy - (itnew->energy)*(itnew->ry);
          Fz = Fz - (itnew->energy)*(itnew->rz);
        }
      }
    }
  }
  // delete irradiator;
  if(savetofile==0)
  {
    printf("%f %f %f\n",Fx,Fy,Fz);
  }
  else
  {
    printf("%f,%f,%f,%f,%f,%f,%f,%f\n",theta,beta,Fx,Fy,Fz,dir_x,dir_y,dir_z);
    fprintf(out,"%f,%f,%f,%f,%f,%f,%f,%f\n",theta,beta,Fx,Fy,Fz,dir_x,dir_y,dir_z);
  }
}





int main(int argc, char* argv[])
{
  #if TEST==1
  cout << "Single triangle single ray tests\n";
  test_1();
  test_2();
  test_3();
  test_4();
  test_5();
  cout << "Multiple triangle single ray tests\n";
  test_6();
  cout << "List tests\n";
  test_7();
  cout << "Shape tests\n";
  test_8();
  test_9();
  cout << "Irradiator tests\n";
  test_10();
  cout << "Raytracing tests\n";
  test_11();
  test_12();
  test_13();
  test_14();
  #endif
  if(strcmp(argv[1],"help")==0)
  {
    printf("Parameters: geometry_file_name direction_x direction_y direction_z distance radius count\n");
    return 0;
  }
  //trace(char * input,double dir_x, double dir_y, double dir_z,double distance,double radius,double count)
  //trace(argv[1],atof(argv[2]), atof(argv[3]), atof(argv[4]),atof(argv[5]),atof(argv[6]),atoi(argv[7]),0);

  int Ntheta = 10;
  int Nbeta = 10;
  int Nmore = 1;
  double theta = 0.0;
  double beta = 0.0;

  Shape * shape = new Shape();
  shape->Load(argv[1]);

  FILE *out = fopen( "output.csv", "w" );
  fprintf(out,"\"theta\",\"beta\",\"FX\",\"FY\",\"FZ\",\"SX\",\"SY\",\"SZ\"\n");
  for(int i = 0; i<(2*Ntheta+1); ++i)
    for(int j = 1; j<Nbeta; ++j)
    {
      theta = M_PI*i/Ntheta;
      beta =  M_PI*j/(2*Nbeta);
      for(int k=0; k<Nmore; ++k)
      {
        trace2(shape,theta,beta,atof(argv[5]),atof(argv[6]),atoi(argv[7]),1,out); 
      }
    }
  fclose(out);
  return 0;
}

