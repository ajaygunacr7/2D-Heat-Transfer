#include <iostream>
#include <vector>
#include <fstream>
#include <string.h>
#include <math.h>

int Problem = 1;

#include "Class.h"
#include "Functions.H"
#include "Matrix.h"
#include "GaussSeidel.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main()
{

  double dR, dPhi;

  int N_phi = 2;
  int N_R = 2;

  std::cout << "Number Of Divisions in Phi direction: ";
  std::cin >> N_phi ;
  std::cout << "\nNumber Of Divisions in Radial direction: ";
  std::cin >> N_R;
  std::cout << "\nSelect Problem(1/2): ";
  std::cin >> Problem;


  dR = (double)   5/(N_R -1)  ;
  dPhi = (double) 360/(N_phi-1);

  std::ofstream OutFile;
  OutFile.open("Dets.txt");
  OutFile << Problem << std::endl << N_phi << std::endl << N_R;

  dR /= 100;
  Mesh.GenerateMesh(dR, dPhi);

  std::vector<std::vector<double>> CoeffMatrix(Mesh.Point.size(),std::vector<double>(Mesh.Point.size(),0.00d));

  std::vector<double> SolutionMatrix(Mesh.Point.size(),400.0d + 273.0d);

  float DeltaT = 10.0f;

  long double AlphaT = 1.16e-5 * DeltaT;
  long double Coeff  = (long double) 2 * dR * 50 / 45;
  long double CoeffWater  = (long double) 2 * dR * 250 / 45;

  if (Problem == 1)
  {
    CoeffWater = 0.0d;
  }

  CreateCoeffMatrix(CoeffMatrix,Mesh,AlphaT,Coeff,CoeffWater);

  double Time = 0;
  unsigned int Int = 0 ;

  double MaxTemp = 201.0d + 273.0d;

  DisplayCoeffMatrix(CoeffMatrix,SolutionMatrix);

  while(MaxTemp >= 200.0d + 273.0d)
  {
      MaxTemp = Mesh.Point[0].Temperature;
      Int++;
      Mesh.CopyMesh(SolutionMatrix);

      ApplyBC(SolutionMatrix,Mesh,Coeff*AlphaT,CoeffWater*AlphaT);

      GaussSeidel(CoeffMatrix,SolutionMatrix,Mesh);

      Time += DeltaT;

    if (Int % 100 == 0)
      std::cout << "Time : " << Time << " Temp : " << MaxTemp <<std::endl;

    for (unsigned int i = 0; i < Mesh.Point.size();++i)
    {
      if (Mesh.Point[i].Temperature > MaxTemp)
        MaxTemp = Mesh.Point[i].Temperature;
    }
  }
  std::cout << "Time Taken " << Time << std::endl;

  WriteData(Mesh);

  std::ofstream File;
  if (Problem == 1)
  {
      File.open("Data1",std::ios_base::app);
  }
  else
  {
      File.open("Data2",std::ios_base::app);
  }
  File << N_phi << '\t' << N_R << '\t' << Time << std::endl;


}
