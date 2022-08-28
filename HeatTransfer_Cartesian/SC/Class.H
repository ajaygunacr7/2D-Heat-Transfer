#pragma once


/*---------------------------------------------------------------------*\
                                STRUCT
\*---------------------------------------------------------------------*/

struct Point
{
  double x;
  double y;
  int z;                      // 2D case Z is set to 0;
  double Temperature;
};

/*---------------------------------------------------------------------*\
                                CLASS
\*---------------------------------------------------------------------*/

class MeshClass
{
public:
  double Delta;
  std::vector<Point> Points;
  long int Length;

public:
  MeshClass(const double);

};

/*---------------------------------------------------------------------*\
                        METHOD DECLARATAIONS
\*---------------------------------------------------------------------*/

MeshClass::MeshClass(const double Del)
{

  Length = 0;

  Delta = Del;

  for (double i = 0.00000d; i <  1.000001d; i+=Del)
  {
     // Creating and Initializing Mesh Data with different points and Temperature = 273K

    for(double j = 0.000000d; j < 1.000001d; j+=Del)
      Points.push_back({j,i,0,0.0d});

    // Number of points in 1 direction (Symmertry Case)
    Length++;
  }

}
/*=======================================================================*/
