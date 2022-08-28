#ifndef GAUSSSEIDEL_H
# define GAUSSSEIDEL_H
/*---------------------------------------------------------------------*\
                          GAUSS SEIDEL SOLVER
\*---------------------------------------------------------------------*/
template<typename T,typename Y, typename X>
void GaussSeidel(T& Matrix, Y& SolutionMatrix, X& Mesh)
{
    long int TotalLength = SolutionMatrix.size();
    long int Iteration = 0;
    double MaxResidue = 0, Residue = 0;
    double Temp = 0.0d;
    double SOR = 1.5d;
    while(true)
      {
        MaxResidue = 0.0000d;
        Iteration++;
        if ( Problem == 0)
        {
          for (long int i = 0; i < TotalLength; ++i)
          {
            double Temperature = SolutionMatrix[i];
            Temp = Mesh.Point[i].Temperature;

            if(i == 0)
            {
              Temperature -= Mesh.Point[i+1].Temperature*Matrix[i][i+1];
            }
            // Outer Points
            else if ( i == Mesh.Nr)
            {
              Temperature -= Mesh.Point[i-1].Temperature*Matrix[i][i-1];
            }
            // Internal Points
            else
            {
              Temperature -= Mesh.Point[i-1].Temperature*Matrix[i][i-1];
              Temperature -= Mesh.Point[i+1].Temperature*Matrix[i][i+1];
            }

            Mesh.Point[i].Temperature = (double) Temperature/Matrix[i][i];
            // Over Relaxation
            Mesh.Point[i].Temperature = Temp + SOR * (Mesh.Point[i].Temperature - Temp);


            // Residue Calculation
            Residue = abs(Temp - Mesh.Point[i].Temperature);
            if(Residue > MaxResidue)
              MaxResidue = Residue;
          }
        }


        else
        {
          for (long int i = 0; i < TotalLength; ++i)
          {
            double Temperature = SolutionMatrix[i];
            Temp = Mesh.Point[i].Temperature;
            for (long int j = 0; j < TotalLength; ++j)
            {
              if(j != i)
              {
                Temperature -= Mesh.Point[j].Temperature*Matrix[i][j];
              }

            }

            Mesh.Point[i].Temperature = (double) Temperature/Matrix[i][i];
            // Over Relaxation
            Mesh.Point[i].Temperature = Temp + SOR * (Mesh.Point[i].Temperature - Temp);

            // Residue Calculation
            Residue = abs(Temp - Mesh.Point[i].Temperature);
            if(Residue > MaxResidue)
              MaxResidue = Residue;
          }

        }


        //if(Iteration%1000 == 0)
        //{
        //  std::cout << "Iteration :" << Iteration << ".\t Residue : " << MaxResidue << std::endl;
        //  WriteData(Mesh,Iteration);
        //}
        //std::cout << Iteration << std::endl;


        if(MaxResidue < 0.00001d)
        {
          //std::cout << Mesh.Point[0].Temperature;
          //WriteData(Mesh,Iteration);
          //std::cout << "\nCovergence Achieved. Total iterations taken : " << Iteration <<" With SOR = " << SOR << std::endl;
          break;

        }

      }
}

#endif

/*---------------------------------------------------------------------*\
                        OPTIMIZED (FAILED FOR NOW)
\*---------------------------------------------------------------------*/



            /*double Temperature = SolutionMatrix[i];
            Temp = Mesh.Point[i].Temperature;

            if(i == 0) //Near the end where lenght 1+length fails
            {
              //Temperature -= Mesh.Point[i+1].Temperature*Matrix[i][i+1];
              //Temperature -= Mesh.Point[i-1].Temperature*Matrix[i][i-1];
              //Temperature -= Mesh.Point[i+Mesh.Nr+1].Temperature*Matrix[i][i+Mesh.Nr+1];
              //Temperature -= Mesh.Point[Mesh.Point.size() - Mesh.Nr].Temperature*Matrix[i][Mesh.Point.size() - Mesh.Nr];
              for(unsigned int j = 1; j < Mesh.Point.size(); j += Mesh.Nr)
              {

                //Temperature -= Mesh.Point[j].Temperature*Matrix[i][j];
              }
            }
            //First Line
            else if ( i > 0 && i < Mesh.Nr)
            {
              if(i == 1)
              {
                Temperature -= Mesh.Point[0].Temperature*Matrix[i][0];
              }
              else
              {
                Temperature -= Mesh.Point[i-1].Temperature*Matrix[i][i-1];
              }
              Temperature -= Mesh.Point[i+1].Temperature*Matrix[i][i+1];
              Temperature -= Mesh.Point[i+Mesh.Nr].Temperature*Matrix[i][i+Mesh.Nr];
              Temperature -= Mesh.Point[TotalLength - Mesh.Nr + i%Mesh.Nr].Temperature*Matrix[i][TotalLength - Mesh.Nr + i%Mesh.Nr];

            }
            //Outer Points
            else if ( i % Mesh.Nr == 0)
            {
              Temperature -= Mesh.Point[i-1].Temperature*Matrix[i][i-1];

              if (i == Mesh.Nr)
              {
                Temperature -= Mesh.Point[i+Mesh.Nr].Temperature*Matrix[i][i+Mesh.Nr];
                Temperature -= Mesh.Point[TotalLength-1].Temperature*Matrix[i][TotalLength-1];
              }
              else if (i == TotalLength-1)
              {
                Temperature -= Mesh.Point[Mesh.Nr].Temperature*Matrix[i][Mesh.Nr];
                Temperature -= Mesh.Point[i-Mesh.Nr].Temperature*Matrix[i][i-Mesh.Nr];
              }
              else
              {
                Temperature -= Mesh.Point[i+Mesh.Nr].Temperature*Matrix[i][i+Mesh.Nr];
                Temperature -= Mesh.Point[i-Mesh.Nr].Temperature*Matrix[i][i-Mesh.Nr];
              }
            }
            //Last Line
            else if( i > TotalLength - Mesh.Nr && i != TotalLength-1)
            {
              if(i%Mesh.Nr == 1)
              {
                Temperature -= Mesh.Point[0].Temperature*Matrix[i][0];
              }
              else
              {
                Temperature -= Mesh.Point[i-1].Temperature*Matrix[i][i-1];
              }
              Temperature -= Mesh.Point[i+1].Temperature*Matrix[i][i+1];
              Temperature -= Mesh.Point[i+Mesh.Nr].Temperature*Matrix[i][i+Mesh.Nr];
              Temperature -= Mesh.Point[i % Mesh.Nr].Temperature*Matrix[i][i%Mesh.Nr];
            }
            //Inner Circle
            else if(i > Mesh.Nr && i < TotalLength - Mesh.Nr && i%Mesh.Nr == 1)
            {
              Temperature -= Mesh.Point[0].Temperature*Matrix[i][0];
              Temperature -= Mesh.Point[i+1].Temperature*Matrix[i][i+1];
              Temperature -= Mesh.Point[i+Mesh.Nr].Temperature*Matrix[i][i+Mesh.Nr];
              Temperature -= Mesh.Point[i-Mesh.Nr].Temperature*Matrix[i][i-Mesh.Nr];
            }
            //Internal Points
            else
            {
              Temperature -= Mesh.Point[i-1].Temperature*Matrix[i][i-1];
              Temperature -= Mesh.Point[i+1].Temperature*Matrix[i][i+1];
              Temperature -= Mesh.Point[i+Mesh.Nr].Temperature*Matrix[i][i+Mesh.Nr];
              Temperature -= Mesh.Point[i-Mesh.Nr].Temperature*Matrix[i][i-Mesh.Nr];
            }

            Mesh.Point[i].Temperature = (double) Temperature/Matrix[i][i];
            // Over Relaxation
            Mesh.Point[i].Temperature = Temp + SOR * (Mesh.Point[i].Temperature - Temp);

            //double


            //std::cout << N << std::endl;


            //Mesh.Point[0].Temperature = Mesh.Point[1].Temperature;
            // Residue Calculation
            Residue = abs(Temp - Mesh.Point[i].Temperature);
            if(Residue > MaxResidue)
              MaxResidue = Residue;
    */
              //Temperature = 0; //Mesh.Point[0].Temperature;
              //int N = 0;
              //for(unsigned int j = 1; j < Mesh.Point.size(); j += Mesh.Nr)
              //{
              //  N++;
              //  Temperature += Mesh.Point[j].Temperature;
              //}
              //Mesh.Point[0].Temperature = (double) Temperature/N;
