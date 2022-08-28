#ifndef GAUSSSEIDEL_H
# define GAUSSSEIDEL_H
/*---------------------------------------------------------------------*\
                          GAUSS SEIDEL SOLVER
\*---------------------------------------------------------------------*/
template<typename T,typename Y, typename X>
void GaussSeidel(const T& Matrix, const Y& SolutionMatrix, X& Mesh, double SOR = 1.0d)
{
    long int TotalLength = Matrix.size();
    long int Iteration = 0;
    double MaxResidue = 0, Residue = 0;
    double Temp = 0.0d;
    //double SOR = 1.957d;
    while(true)
      {
        MaxResidue = 0.0000d;
        Iteration++;
        for (long int i = Mesh.Length; i < TotalLength; ++i)
        {
          double Temperature = SolutionMatrix[i];
          Temp = Mesh.Points[i].Temperature;

          if(TotalLength - i < Mesh.Length && i != TotalLength - 1) //Near the end where lenght 1+length fails
          {
            Temperature -= Mesh.Points[i+1].Temperature*Matrix[i][i+1];
            Temperature -= Mesh.Points[i-1].Temperature*Matrix[i][i+1];
            Temperature -= Mesh.Points[i-Mesh.Length].Temperature*Matrix[i][i-Mesh.Length];
            Mesh.Points[i].Temperature = (double)Temperature/Matrix[i][i];
          }
          else if (i == TotalLength - 1) // last Point
          {
            Temperature -= Mesh.Points[i-1].Temperature*Matrix[i][i-1];
            Temperature -= Mesh.Points[i-Mesh.Length].Temperature*Matrix[i][i-Mesh.Length];
            Mesh.Points[i].Temperature = (double)Temperature/Matrix[i][i];
          }
          else
          {
            Temperature -= Mesh.Points[i+Mesh.Length].Temperature*Matrix[i][i+Mesh.Length];
            Temperature -= Mesh.Points[i+1].Temperature*Matrix[i][i+1];
            Temperature -= Mesh.Points[i-1].Temperature*Matrix[i][i-1];
            Temperature -= Mesh.Points[i-Mesh.Length].Temperature*Matrix[i][i-Mesh.Length];
            Mesh.Points[i].Temperature = (double)Temperature/Matrix[i][i];
          }

          // Over Relaxation
          Mesh.Points[i].Temperature = Temp + SOR * (Mesh.Points[i].Temperature - Temp);


          // Residue Calculation
          Residue = abs(Temp - Mesh.Points[i].Temperature);
          if(Residue > MaxResidue)
            MaxResidue = Residue;
        }

        if(Iteration%1000 == 0)
        {
          std::cout << "Iteration :" << Iteration << ".\t Residue : " << MaxResidue << std::endl;
          WriteData(Mesh,Iteration);
        }
        if(MaxResidue < 0.00001d)
        {
          WriteData(Mesh,Iteration);
          std::cout << "\nCovergence Achieved. Total iterations taken : " << Iteration <<" With SOR = " << SOR << std::endl;
          break;
        }

      }
}

#endif
