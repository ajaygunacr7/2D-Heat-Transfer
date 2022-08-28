#pragma once
/*---------------------------------------------------------------------*\
                        MATRIX VALUE ASSIGNMENT
\*---------------------------------------------------------------------*/

template<typename T,typename Y>
void CreateCoeffMatrix(T& Matrix,Y& Mesh,double Coeff = 0.0d)
{
  long int Length = Mesh.Length;

  for(unsigned int i = 0; i < Mesh.Points.size(); ++i)
  {
    // North (South) = i +(-) Length, East(West) = i +(-) 1

    //Single Out Right Boundary
     if(i%Length == Length - 1 && i!= Length-1 && i/Length != Length - 1 )
     {

       Matrix[i][i] = 4.0d + Coeff;

       // West
       Matrix[i][i-1] = -2.0d;

       // North, South
       Matrix[i][i+Length] = -1.0d;
       Matrix[i][i-Length] = -1.0d;
     }

     //Single Out Left Boundary
     else if(i % Length == 0 && i!=0 && i/Length != Length - 1 )
     {

       Matrix[i][i] = 4.0d;

       // East
       Matrix[i][i+1] = -2.0d;

       // North, South
       Matrix[i][i+Length] = -1.0d;
       Matrix[i][i-Length] = -1.0d;
     }

     //Single Out Top Boundary
     else if(i/Length == Length - 1 && i % Length != 0 && i%Length != Length - 1)
     {
       Matrix[i][i] = 4.0d;
       // East, West
       Matrix[i][i+1] = -1.0d;
       Matrix[i][i-1] = -1.0d;

       // South
       Matrix[i][i-Length] = -2.0d;
     }

     //Single Out Bottom Boundary
     else if(i / Length == 0 && i != 0 && i%Length != Length - 1)
     {
       Matrix[i][i] = 1.0d;
       //// East, West
       //Matrix[i][i+1] = -1;
       //Matrix[i][i-1] = -1;

       //// North
       //Matrix[i][i+Length] = -2;
     }

     //Bottom Left Corner
     else if (i == 0)
     {
       Matrix[i][i] = 1.0d;

       //East
       //Matrix[i][i+1] = -2;

       ////North
       //Matrix[i][i+Length] = -2;
     }

     // Bottom Right Corner
     else if(i == Length - 1)
     {
       Matrix[i][i] = 1.0d;

       // West
       //Matrix[i][i-1] = -2;

       ////North
       //Matrix[i][i+Length] = -2;
     }

     // Top Right Corner
     else if(i%Length == Length - 1 && i/Length == Length - 1)
     {

       Matrix[i][i] = 4.0d + Coeff;

       // West
       Matrix[i][i-1] = -2.0d;

       //North
       Matrix[i][i-Length] = -2.0d;
     }

     //Top Left Corner
     else if(i % Length == 0 && i/Length == Length - 1)
     {
       Matrix[i][i] = 4.0d;

       // East
       Matrix[i][i+1] = -2.0d;

       //South
       Matrix[i][i-Length] = -2.0d;
     }

     // Points in the Centre
     else
     {
       Matrix[i][i] = 4.0d;

       // East, West
       Matrix[i][i+1] = -1.0d;
       Matrix[i][i-1] = -1.0d;

       // North, South
       Matrix[i][i+Length] = -1.0d;
       Matrix[i][i-Length] = -1.0d;
     }

  }
}

/*---------------------------------------------------------------------*\
                      DISPLAY MATRIX FUNCTION
\*---------------------------------------------------------------------*/


template<typename T,typename Y>
void DisplayCoeffMatrix(const T& Matrix,const Y& SolutionMatrix)
{
  for(unsigned int i = 0; i < Matrix.size(); ++i)
  {
    std::cout << "|";
    for(unsigned int j = 0; j <  Matrix.size() ; ++j)
    {
      if(Matrix[i][j] < 0)
        std::cout << Matrix[i][j] << " ";
      else
        std::cout << " " << Matrix[i][j] << " ";
    }

    std::cout << "| T[" << i+1 << "] = " << SolutionMatrix[i] << std::endl;
  }
}
