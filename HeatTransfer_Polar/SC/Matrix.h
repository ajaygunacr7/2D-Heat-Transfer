#ifndef MATRIX_H
# define MATRIX_H
/*---------------------------------------------------------------------*\
                      CREATE COEFF MATRIX FUNCTION
\*---------------------------------------------------------------------*/

template <typename T, typename Y>
void CreateCoeffMatrix(T& CoeffMatrix, Y& Mesh,long double Coeff = 1.0d, long double HCoeff = 0.0d, long double HCoeffWater = 0.0d)
{
  if (Problem == 0)
  {
    for (unsigned int i = 0; i < Mesh.Point.size(); ++i)
    {
      if (i == 0)
      {
        CoeffMatrix[i][i] += (double)   4 * Coeff / pow(Mesh.DelR,2);

        CoeffMatrix[i][i+1] += (double) -4*Coeff/pow(Mesh.DelR,2);
      }
      //First Line
      else if (i < Mesh.Nr && i != 0)
      {

        CoeffMatrix[i][i] += (double) 2/pow(Mesh.DelR,2);
        CoeffMatrix[i][i] *= Coeff;

        CoeffMatrix[i][i-1] += (double) 1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i-1] += (double) -1/pow(Mesh.DelR,2);
        CoeffMatrix[i][i-1] *= Coeff;

        CoeffMatrix[i][i+1] += (double) -1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i+1] += (double) -1/pow(Mesh.DelR,2);
        CoeffMatrix[i][i+1] *= Coeff;

      }

      else // if (i == Mesh.Nr)
      {

        CoeffMatrix[i][i] += (double)  2/pow(Mesh.DelR,2);
        CoeffMatrix[i][i] += (double) HCoeff/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i] += (double) HCoeff/pow(Mesh.DelR,2);
        CoeffMatrix[i][i] *= Coeff;

        //CoeffMatrix[i][i-1] += (double) 1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i-1] += (double) -1/pow(Mesh.DelR,2);
        //CoeffMatrix[i][i-1] += (double) -1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i-1] += (double) -1/pow(Mesh.DelR,2);
        CoeffMatrix[i][i-1] *= Coeff;

      }
      CoeffMatrix[i][i] += 1;
    }
  }

  ///////////////////////////
  //  2d Problem of A
  ///////////////////////////


  else if (Problem == 1)
  {
    for (unsigned int i = 0; i < Mesh.Point.size(); ++i)
    {
      if (i == 0)
      {
          for(unsigned int j = 1; j < Mesh.Point.size(); j += Mesh.Nr)
          {
          CoeffMatrix[i][i] += (double)   4 * Coeff / pow(Mesh.DelR,2);


          CoeffMatrix[i][j] += (double) -4*Coeff/pow(Mesh.DelR,2);

          //CoeffMatrix[i][Mesh.Point.size() - Mesh.Nr] *= Coeff;

          //CoeffMatrix[i][1] = CoeffMatrix[i][Mesh.Point.size() - Mesh.Nr];
        }
        CoeffMatrix[i][i] += (360/90 - 1);
      }
      //First Line
      else if (i < Mesh.Nr && i != 0)
      {

        CoeffMatrix[i][i] += (double) 2/pow(Mesh.DelR,2);
        CoeffMatrix[i][i] += (double) 2/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
        CoeffMatrix[i][i] *= Coeff;

        CoeffMatrix[i][i-1] += (double) 1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i-1] += (double) -1/pow(Mesh.DelR,2);
        CoeffMatrix[i][i-1] *= Coeff;

        CoeffMatrix[i][i+1] += (double) -1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i+1] += (double) -1/pow(Mesh.DelR,2);
        CoeffMatrix[i][i+1] *= Coeff;


        CoeffMatrix[i][i + Mesh.Nr] += (double) -1/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
        CoeffMatrix[i][i + Mesh.Nr] *= Coeff;

        unsigned int Pos = Mesh.Point.size() - Mesh.Nr + (i % Mesh.Nr);

        CoeffMatrix[i][Pos] = CoeffMatrix[i][i + Mesh.Nr];
      }



      //Outer Circle
      else if ( i % Mesh.Nr == 0)
      {

        CoeffMatrix[i][i] += (double)  2/pow(Mesh.DelR,2);
        CoeffMatrix[i][i] += (double) 2/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));

        double HTCoeff = HCoeff;
        if (HCoeffWater != 0.0d)
        {
          if(Mesh.Point[i].Phi < 3.14059268)  // delibrately made it 0
            HTCoeff = HCoeff;
          else
            HTCoeff = HCoeffWater;
        }


        CoeffMatrix[i][i] += (double) HTCoeff/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i] += (double) HTCoeff/pow(Mesh.DelR,2);

        CoeffMatrix[i][i] *= Coeff;

        CoeffMatrix[i][i-1] += (double) 1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i-1] += (double) -1/pow(Mesh.DelR,2);

        CoeffMatrix[i][i-1] += (double) -1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i-1] += (double) -1/pow(Mesh.DelR,2);

        CoeffMatrix[i][i-1] *= Coeff;

        if(i == Mesh.Nr)
        {
          CoeffMatrix[i][i+Mesh.Nr] += (double) -1/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
          CoeffMatrix[i][i+Mesh.Nr] *= Coeff;

          unsigned int Pos = Mesh.Point.size()-Mesh.Nr + (i % Mesh.Nr);

          CoeffMatrix[i][Pos] = CoeffMatrix[i][i+Mesh.Nr];
        }

        else if(i == Mesh.Point.size() - 1)
        {
          CoeffMatrix[i][i-Mesh.Nr] += (double) -1/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
          CoeffMatrix[i][i-Mesh.Nr] *= Coeff;

          unsigned int Pos = (i % Mesh.Nr);

          CoeffMatrix[i][Pos] = CoeffMatrix[i][i-Mesh.Nr];
        }
        else
        {
          CoeffMatrix[i][i-Mesh.Nr] += (double) -1/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
          CoeffMatrix[i][i-Mesh.Nr] *= Coeff;
          CoeffMatrix[i][i+Mesh.Nr] = CoeffMatrix[i][i-Mesh.Nr];
        }
      }
      //Last Line
      else if (i > Mesh.Point.size() - Mesh.Nr && i != Mesh.Point.size() - 1)
      {

        CoeffMatrix[i][i] += (double)  2/pow(Mesh.DelR,2);
        CoeffMatrix[i][i] += (double) 2/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
        CoeffMatrix[i][i] *= Coeff;

        CoeffMatrix[i][i-1] += (double) 1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i-1] += (double) -1/pow(Mesh.DelR,2);
        CoeffMatrix[i][i-1] *= Coeff;

        CoeffMatrix[i][i+1] += (double) -1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i+1] += (double) -1/pow(Mesh.DelR,2);
        CoeffMatrix[i][i+1] *= Coeff;

        CoeffMatrix[i][i-Mesh.Nr] = (double) -1/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
        CoeffMatrix[i][i-Mesh.Nr] *= Coeff;

        unsigned int Pos = (i % Mesh.Nr);

        CoeffMatrix[i][Pos] = CoeffMatrix[i][i-Mesh.Nr];
      }

      // Inner Circle
      else if (i % Mesh.Nr == 1 && i != 1)
      {

        CoeffMatrix[i][i] += (double) 2/pow(Mesh.DelR,2);
        CoeffMatrix[i][i] += (double) 2/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
        CoeffMatrix[i][i] *= Coeff;

        CoeffMatrix[i][0] += (double) 1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][0] += (double) -1/pow(Mesh.DelR,2);
        CoeffMatrix[i][0] *= Coeff;

        CoeffMatrix[i][i+1] += (double) -1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i+1] += (double) -1/pow(Mesh.DelR,2);
        CoeffMatrix[i][i+1] *= Coeff;



        if(i == Mesh.Point.size() - Mesh.Nr)
        {

          CoeffMatrix[i][i-Mesh.Nr]  += (double) -1/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
          CoeffMatrix[i][i-Mesh.Nr] *= Coeff;

          unsigned int Pos = 1;

          CoeffMatrix[i][Pos] = CoeffMatrix[i][i-Mesh.Nr];
        }
        else
        {
          CoeffMatrix[i][i+Mesh.Nr] = (double) -1/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
          CoeffMatrix[i][i+Mesh.Nr] *= Coeff;

          CoeffMatrix[i][i-Mesh.Nr] = CoeffMatrix[i][i+Mesh.Nr];
        }

      }
      else
      {

        CoeffMatrix[i][i] +=  (double) 2/pow(Mesh.DelR,2);
        CoeffMatrix[i][i] += (double) 2/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
        CoeffMatrix[i][i] *= Coeff;

        CoeffMatrix[i][i-1] +=  (double) 1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i-1] += (double) -1/pow(Mesh.DelR,2);
        CoeffMatrix[i][i-1] *= Coeff;

        CoeffMatrix[i][i+1] +=  (double) -1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i+1] += (double) -1/pow(Mesh.DelR,2);
        CoeffMatrix[i][i+1] *= Coeff;

        CoeffMatrix[i][i+Mesh.Nr] = (double) -1/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
        CoeffMatrix[i][i+Mesh.Nr] *= Coeff;

        CoeffMatrix[i][i-Mesh.Nr] = CoeffMatrix[i][i+Mesh.Nr];
      }
      CoeffMatrix[i][i] += 1;
    }
  }

  ///////////////////////////
  //  2d Problem of B
  ///////////////////////////
  else
  {
    for (unsigned int i = 0; i < Mesh.Point.size(); ++i)
    {
      if (i == 0)
      {
          for(unsigned int j = 1; j < Mesh.Point.size(); j += Mesh.Nr)
          {
          CoeffMatrix[i][i] += 1;


          CoeffMatrix[i][j] -= 1;

          //CoeffMatrix[i][Mesh.Point.size() - Mesh.Nr] *= Coeff;

          //CoeffMatrix[i][1] = CoeffMatrix[i][Mesh.Point.size() - Mesh.Nr];
        }
        CoeffMatrix[i][i] -= 1;
      }
      //First Line
      else if (i < Mesh.Nr && i != 0)
      {

        CoeffMatrix[i][i] += (double) 2/pow(Mesh.DelR,2);
        CoeffMatrix[i][i] += (double) 2/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
        CoeffMatrix[i][i] *= Coeff;

        CoeffMatrix[i][i-1] += (double) 1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i-1] += (double) -1/pow(Mesh.DelR,2);
        CoeffMatrix[i][i-1] *= Coeff;

        CoeffMatrix[i][i+1] += (double) -1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i+1] += (double) -1/pow(Mesh.DelR,2);
        CoeffMatrix[i][i+1] *= Coeff;


        CoeffMatrix[i][i + Mesh.Nr] += (double) -1/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
        CoeffMatrix[i][i + Mesh.Nr] *= Coeff;

        unsigned int Pos = Mesh.Point.size() - Mesh.Nr + (i % Mesh.Nr);

        CoeffMatrix[i][Pos] = CoeffMatrix[i][i + Mesh.Nr];
      }



      //Outer Circle
      else if ( i % Mesh.Nr == 0)
      {

        CoeffMatrix[i][i] += (double)  2/pow(Mesh.DelR,2);
        CoeffMatrix[i][i] += (double) 2/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));

        double HTCoeff;
        if (HCoeffWater != 0.0d)
        {
          if(Mesh.Point[i].Phi < 3.14059268)  // delibrately made it 0
            HTCoeff = HCoeff;
          else
            HTCoeff = HCoeffWater;
        }
        else
          HTCoeff = HCoeff;

        CoeffMatrix[i][i] += (double) HTCoeff/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i] += (double) HTCoeff/pow(Mesh.DelR,2);

        CoeffMatrix[i][i] *= Coeff;

        CoeffMatrix[i][i-1] += (double) 1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i-1] += (double) -1/pow(Mesh.DelR,2);

        CoeffMatrix[i][i-1] += (double) -1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i-1] += (double) -1/pow(Mesh.DelR,2);

        CoeffMatrix[i][i-1] *= Coeff;

        if(i == Mesh.Nr)
        {
          CoeffMatrix[i][i+Mesh.Nr] += (double) -1/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
          CoeffMatrix[i][i+Mesh.Nr] *= Coeff;

          unsigned int Pos = Mesh.Point.size()-Mesh.Nr + (i % Mesh.Nr);

          CoeffMatrix[i][Pos] = CoeffMatrix[i][i+Mesh.Nr];
        }

        else if(i == Mesh.Point.size() - 1)
        {
          CoeffMatrix[i][i-Mesh.Nr] += (double) -1/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
          CoeffMatrix[i][i-Mesh.Nr] *= Coeff;

          unsigned int Pos = (i % Mesh.Nr);

          CoeffMatrix[i][Pos] = CoeffMatrix[i][i-Mesh.Nr];
        }
        else
        {
          CoeffMatrix[i][i-Mesh.Nr] += (double) -1/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
          CoeffMatrix[i][i-Mesh.Nr] *= Coeff;
          CoeffMatrix[i][i+Mesh.Nr] = CoeffMatrix[i][i-Mesh.Nr];
        }
      }
      //Last Line
      else if (i > Mesh.Point.size() - Mesh.Nr && i != Mesh.Point.size() - 1)
      {

        CoeffMatrix[i][i] += (double)  2/pow(Mesh.DelR,2);
        CoeffMatrix[i][i] += (double) 2/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
        CoeffMatrix[i][i] *= Coeff;

        CoeffMatrix[i][i-1] += (double) 1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i-1] += (double) -1/pow(Mesh.DelR,2);
        CoeffMatrix[i][i-1] *= Coeff;

        CoeffMatrix[i][i+1] += (double) -1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i+1] += (double) -1/pow(Mesh.DelR,2);
        CoeffMatrix[i][i+1] *= Coeff;

        CoeffMatrix[i][i-Mesh.Nr] = (double) -1/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
        CoeffMatrix[i][i-Mesh.Nr] *= Coeff;

        unsigned int Pos = (i % Mesh.Nr);

        CoeffMatrix[i][Pos] = CoeffMatrix[i][i-Mesh.Nr];
      }

      // Inner Circle
      else if (i % Mesh.Nr == 1 && i != 1)
      {

        CoeffMatrix[i][i] += (double) 2/pow(Mesh.DelR,2);
        CoeffMatrix[i][i] += (double) 2/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
        CoeffMatrix[i][i] *= Coeff;

        CoeffMatrix[i][0] += (double) 1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][0] += (double) -1/pow(Mesh.DelR,2);
        CoeffMatrix[i][0] *= Coeff;

        CoeffMatrix[i][i+1] += (double) -1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i+1] += (double) -1/pow(Mesh.DelR,2);
        CoeffMatrix[i][i+1] *= Coeff;



        if(i == Mesh.Point.size() - Mesh.Nr)
        {

          CoeffMatrix[i][i-Mesh.Nr]  += (double) -1/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
          CoeffMatrix[i][i-Mesh.Nr] *= Coeff;

          unsigned int Pos = 1;

          CoeffMatrix[i][Pos] = CoeffMatrix[i][i-Mesh.Nr];
        }
        else
        {
          CoeffMatrix[i][i+Mesh.Nr] = (double) -1/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
          CoeffMatrix[i][i+Mesh.Nr] *= Coeff;

          CoeffMatrix[i][i-Mesh.Nr] = CoeffMatrix[i][i+Mesh.Nr];
        }

      }
      else
      {

        CoeffMatrix[i][i] +=  (double) 2/pow(Mesh.DelR,2);
        CoeffMatrix[i][i] += (double) 2/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
        CoeffMatrix[i][i] *= Coeff;

        CoeffMatrix[i][i-1] +=  (double) 1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i-1] += (double) -1/pow(Mesh.DelR,2);
        CoeffMatrix[i][i-1] *= Coeff;

        CoeffMatrix[i][i+1] +=  (double) -1/(Mesh.Point[i].R * 2 * Mesh.DelR);
        CoeffMatrix[i][i+1] += (double) -1/pow(Mesh.DelR,2);
        CoeffMatrix[i][i+1] *= Coeff;

        CoeffMatrix[i][i+Mesh.Nr] = (double) -1/(pow(Mesh.Point[i].R,2)*pow(Mesh.DelPhi,2));
        CoeffMatrix[i][i+Mesh.Nr] *= Coeff;

        CoeffMatrix[i][i-Mesh.Nr] = CoeffMatrix[i][i+Mesh.Nr];
      }
      CoeffMatrix[i][i] += 1;
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
#endif
