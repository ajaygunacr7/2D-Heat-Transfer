#pragma once

/*---------------------------------------------------------------------*\
                        FUNCTION DEFINITIONS
\*---------------------------------------------------------------------*/

template<typename T, typename Y>
void ApplyBC(T& SolutionMatrix,Y& Mesh, long double Coeff = 0.0d, long double CoeffW = 0.0d)
{
  if (Problem == 1)
  {
    SolutionMatrix[0] *= 360/90;
    for(unsigned int i = Mesh.Nr; i < SolutionMatrix.size(); i += Mesh.Nr)
    {
      //std::cout << SolutionMatrix[i] << std::endl;
      double value  = (double) Coeff*293.0d*(1/(0.05 * 2 * Mesh.DelR) + 1/(pow(Mesh.DelR,2)));
      SolutionMatrix[i] += value;

      //std::cout << value<< std::endl;
    }
  }
  else
  {
    SolutionMatrix[0] = 0;
    for(unsigned int i = Mesh.Nr; i < SolutionMatrix.size(); i += Mesh.Nr)
    {
      //std::cout << SolutionMatrix[i] << std::endl;
      double HTCoeff;
      if (CoeffW != 0.0d)
      {
        if(Mesh.Point[i].Phi < 3.14059268)  // delibrately made it 0
          HTCoeff = Coeff*293.0d;
        else
          HTCoeff = CoeffW*278.0d;
      }
      double value  = (double) HTCoeff*(1/(0.05 * 2 * Mesh.DelR) + 1/(pow(Mesh.DelR,2)));
      SolutionMatrix[i] += value;

      //std::cout << value<< std::endl;
    }
  }
}

/*=======================================================================*/

template<typename T>
void WriteData(T& Mesh)

{
  std::string Name = "Temp.txt";
  std::ofstream File(Name.c_str());

  // Giving Field Names To Coloumns
  //File << "x,y,z,Temperature" << std::endl;

  // Writing Files to data 'X, Y, Z, Temperature' Format
  for(unsigned int i = 0; i < Mesh.Point.size(); ++i)
  {
    double x, y;
    x = Mesh.Point[i].R * cos(Mesh.Point[i].Phi)*100;
    y = Mesh.Point[i].R * sin(Mesh.Point[i].Phi)*100;
    File << Mesh.Point[i].Temperature << std::endl;
  }
}

/*=======================================================================*/
