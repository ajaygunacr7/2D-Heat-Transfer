#pragma once

/*---------------------------------------------------------------------*\
                        FUNCTION DEFINITIONS
\*---------------------------------------------------------------------*/
template<typename T,typename Y>
void ImposeBoundaryCondition(T& Mesh, Y& SolutionMatrix, double Coeff = 0.0d)
{
  // Temperature = 100*sin(pi*x)
  for(unsigned int i = 0; i < Mesh.Length; ++i)
  {
    Mesh.Points[i].Temperature = 0;//(100.0d*sin(3.14159d*Mesh.Points[i].x)) + 273.0d;
    SolutionMatrix[i] = 1.0d;
  }
  if (Coeff != 0.0d)
  {
    for(unsigned int i = 2*(Mesh.Length) - 1; i < SolutionMatrix.size(); i+=Mesh.Length)
      SolutionMatrix[i] = Coeff; //pow(Mesh.Delta,2);
  }
}

/*=======================================================================*/
template<typename T>
void WriteData(T& Mesh,const int Iteration)

{
  std::string Name = "PlotData/mesh" + std::to_string(Iteration) + ".csv";
  std::ofstream File(Name.c_str());

  // Giving Field Names To Coloumns
  File << "x,y,z,Temperature" << std::endl;

  // Writing Files to data 'X, Y, Z, Temperature' Format
  for(unsigned int i = 0; i < Mesh.Points.size(); ++i)
  {
    File << Mesh.Points[i].x << "," << Mesh.Points[i].y << "," << Mesh.Points[i].z << ","  << Mesh.Points[i].Temperature << std::endl;
  }
}

/*=======================================================================*/
