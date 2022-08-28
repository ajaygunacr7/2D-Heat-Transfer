#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <sstream>

/*---------------------------------------------------------------------*\
                            INCLUDES
\*---------------------------------------------------------------------*/

#include "Class.H"          // Contains the Mesh Class and its Methods
#include "Functions.H"      // Contains Functions
#include "Matrix.H"         // Contains Functions Related to Matrix
#include "GaussSeidel.h"    // Contains Gauss Seildel Solver


/*---------------------------------------------------------------------*\
                            MAIN FUNCTION
\*---------------------------------------------------------------------*/

int main()
{
  int Nodes;

  std::cout << "Enter Number Of Divisions in any one Coorinate Axis (NxN) :";
  std::cin >> Nodes;

  // Delta grid spacing (Thinking of converting this to an array for further use)
  double Del = (double) 1/Nodes;


  // Creating Object "Mesh" and calling the constructor
  MeshClass Mesh(Del);

  std::cout << "\nTotal Points in Mesh : " << Mesh.Points.size() << std::endl;

  // Creating and Initializing Coefficent Matrix with 0s
  std::vector<std::vector<double>> CoeffMatrix(Mesh.Points.size(),std::vector<double> (Mesh.Points.size(),0.0d));


  // Creating and Initializing Solution Matrix with 0
  std::vector<double> SolutionMatrix (Mesh.Points.size(), 0.0d);

  int Option;

  std::cout << "\n Which Problem Do You Want To Solve (1/2)  ";
  std::cin >> Option;

  double SOR;

  std::cout << "\n Enter the SOR : ";
  std::cin >> SOR;

  if (Option == 1)
  {
    // Imposing the Bounary Condition T = 100*sin(pi*x)
    ImposeBoundaryCondition(Mesh,SolutionMatrix);


    // Creating the Coefficient Matrix
    CreateCoeffMatrix(CoeffMatrix,Mesh);
  }
  else
  {

    double AmbientTemperature = 30.0d + 273.0d;

    float H = 20.0f, K = 30.0f;

    double Coeff = (double) 2*Del*H/(K);

    ImposeBoundaryCondition(Mesh,SolutionMatrix, Coeff*AmbientTemperature);

    std::cout << "Coefficent : " << Coeff << std::endl;

    CreateCoeffMatrix(CoeffMatrix,Mesh,Coeff);

  }

  // Print CoeffMatrix
  DisplayCoeffMatrix(CoeffMatrix,SolutionMatrix);

  // Writing data into file for post-processing
  //WriteData(Mesh,0);


  // Calling Gauss Seidel Solver
  GaussSeidel(CoeffMatrix,SolutionMatrix,Mesh,SOR);

}
