#pragma once

struct Pnt
{
  double R, Phi;
  double Temperature;
};

class Msh
{
public:
  std::vector<Pnt> Point;
  unsigned int Nr;
  double DelR, DelPhi;
public:

  void GenerateMesh( const double dR, const double dPhi );

  void DisplayMesh();

  void CopyMesh(std::vector<double>& M);

  void copy(std::vector<double>& R,std::vector<double>& T, std::vector<double>& Phi);

}Mesh,MeshOld;


void Msh::GenerateMesh(const double dR, const double dPhi)
{

  double Pi = (double) acos(0)/90;
  double PhiR;

  DelR = dR;
  DelPhi = dPhi * Pi;

  unsigned int L = 0;
  Point.push_back({0,0,400 + 273.0d});

  if (dPhi == 0)
  {
    for (double Radius = dR; Radius <= 0.050001d; Radius += dR)
    {
      Point.push_back({Radius, 0.0d, 400.0d + 273.0d});
      L++;
    }
  }
  else
  {
    for (double PhiD = 0; PhiD < 360; PhiD += dPhi)
    {
      PhiR = PhiD * Pi;
      L = 0;
      for (double Radius = dR; Radius <= 0.050001d; Radius += dR)
      {
        Point.push_back({Radius, PhiR, 400.0d + 273.0d});
        L++;
      }
    }
  }
  //std::cout << L << std::endl;
  Nr = L;
}

void Msh::DisplayMesh()
{
  for (unsigned int i = 0; i < Point.size(); ++i)
  {

    std::cout << "Point [" << i +1 << "] R = "<< Point[i].R;
    std::cout << " Phi = " << Point[i].Phi << " T = "<< Point[i].Temperature;
    std::cout << std::endl;

  }
}

void Msh::CopyMesh(std::vector<double>& M)
{
  for (unsigned int i = 0; i < Point.size(); ++i)
    M[i] = Point[i].Temperature;
}



void Msh::copy(std::vector<double>& R,std::vector<double>& T, std::vector<double>& Phi)
{
  for (unsigned int i = 0; i < Point.size(); ++i)
  {
    R.push_back(Point[i].R);
    T.push_back(Point[i].Temperature);
    Phi.push_back(Point[i].Phi);
  }
}
