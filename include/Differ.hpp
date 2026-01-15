










# include <iostream>
# include <algorithm>





# include "Array.hpp"
# include "Global.hpp"

using namespace std;
using namespace ARRAY;









double Minmod(double R)
{
    double Minmod = 0.0;

    if(R > 0)
        Minmod = min(R,1.0);
    else
        Minmod = 0.0;
    
    return Minmod;
}

double Van_Leer(double R)
{
    double Van_Leer = 0.0;

    Van_Leer = ( R + abs(R)) / (1.0 + abs(R));

    return Van_Leer;
}

double Superbee(double R)
{
    double Superbee = 0.0;

    Superbee = max({0.0,min(2 * R,1.0),min(R,2.0)});

    return Superbee;
}

double Van_Albada(double R)
{
    double Van_Albada = 0.0;

    Van_Albada = (pow(R,2) + R) / (1.0 + pow(R,2));

    return Van_Albada;
}

double Double_Minmod(double R)
{
    double Double_Minmod = 0.0;

    if(R > 0)
        Double_Minmod = min({2 * R,1.0,(1 + R) / 2});
    else
        Double_Minmod = 0.0;

    return Double_Minmod; 
}












void Diff_Initial(int direction, Array<double,1> &xnode, Array<double,1> &ynode, Array<double,1> &znode, Array<double,4> &FLR, Array<double,3> &fi, int bc, int limit)
{
}












void MUSCL_1(int direction, Array<double,1> &xnode, Array<double,1> &ynode, Array<double,1> &znode, Array<double,4> &FLR, Array<double,3> &fi, int bc, int limit)
{
    int ni = xnode.GetSize() - 2 * bc;
    int nj = ynode.GetSize() - 2 * bc;
    int nk = znode.GetSize() - 2 * bc;

    for(int i = bc - 1;i < ni + bc;i ++)
        for(int j = bc - 1;j < nj + bc;j ++)
            for(int k = bc - 1;k < nk + bc;k ++)
            {
                if(direction == 1)
                {
                    FLR(i,j,k,0) = fi(i,j,k);
                    FLR(i,j,k,1) = fi(i + 1,j,k);
                }
                else if (direction == 2)
                {
                    FLR(i,j,k,0) = fi(i,j,k);
                    FLR(i,j,k,1) = fi(i,j + 1,k);
                }
                else if (direction == 3)
                {
                    FLR(i,j,k,0) = fi(i,j,k);
                    FLR(i,j,k,1) = fi(i,j,k + 1);
                }
            }
}












void MUSCL_2(int direction,Array<double,1> &xnode, Array<double,1> &ynode, Array<double,1> &znode, Array<double,4> &FLR,Array<double,3> &fi, int bc, int limit)
{
    int ni = xnode.GetSize() - 2 * bc;
    int nj = ynode.GetSize() - 2 * bc;
    int nk = znode.GetSize() - 2 * bc;

    double kk = 1.0 / 3.0;
    double A0 = 0.0,A1 = 0.0,A2 = 0.0,L = 0.0,R = 0.0;
    double L1 = 0.0,L2 = 0.0,R1 = 0.0,R2 = 0.0;


    for(int i = bc - 1;i < ni + bc;i ++)
        for(int j = bc - 1;j < nj + bc;j ++)
            for(int k = bc - 1;k < nk + bc;k ++)
            {


                if(direction == 1)
                {
                    A0 = fi(i,j,k) - fi(i - 1,j,k);
                    A1 = fi(i + 1,j,k) - fi(i,j,k);
                    A2 = fi(i + 2,j,k) - fi(i + 1,j,k);
                }
                else if(direction == 2)
                {
                    A0 = fi(i,j,k) - fi(i,j - 1,k);
                    A1 = fi(i,j + 1,k) - fi(i,j,k);
                    A2 = fi(i,j + 2,k) - fi(i,j + 1,k);
                }
                else if(direction == 3)
                {
                    A0 = fi(i,j,k) - fi(i,j,k - 1);
                    A1 = fi(i,j,k + 1) - fi(i,j,k);
                    A2 = fi(i,j,k + 2) - fi(i,j,k + 1);
                }

                L = A1 / (A0 + 1e-12);
                R = A1 / (A2 + 1e-12);

                switch(limit)
                {
                case 0:
                    L1 = Minmod(L);
                    L2 = Minmod(1.0 / L);
                    R1 = Minmod(R);
                    R2 = Minmod(1.0 / R);
                    break;
                case 1:
                    L1 = Van_Leer(L);
                    L2 = Van_Leer(1.0 / L);
                    R1 = Van_Leer(R);
                    R2 = Van_Leer(1.0 / R);
                    break;
                case 2:
                    L1 = Van_Albada(L);
                    L2 = Van_Albada(1.0 / L);
                    R1 = Van_Albada(R);
                    R2 = Van_Albada(1.0 / R);
                    break;
                case 3:
                    L1 = Superbee(L);
                    L2 = Superbee(1.0 / L);
                    R1 = Superbee(R);
                    R2 = Superbee(1.0 / R);
                    break;
                case 4:
                    L1 = Double_Minmod(L);
                    L2 = Double_Minmod(1.0 / L);
                    R1 = Double_Minmod(R);
                    R2 = Double_Minmod(1.0 / R);
                    break;
                }

                if(direction == 1)
                {
                    FLR(i,j,k,0) = fi(i,j,k) + 0.25 * ((1.0 - kk) * L1 + (1 + kk) * L2 * L) * A0;
                    FLR(i,j,k,1) = fi(i + 1,j,k) - 0.25 * ((1.0 - kk) * R1 + (1 + kk) * R2 * R) * A2;
                }
                else if(direction == 2)
                {
                    FLR(i,j,k,0) = fi(i,j,k) + 0.25 * ((1.0 - kk) * L1 + (1 + kk) * L2 * L) * A0;
                    FLR(i,j,k,1) = fi(i,j + 1,k) - 0.25 * ((1.0 - kk) * R1 + (1 + kk) * R2 * R) * A2;
                }
                else if(direction == 3)
                {
                    FLR(i,j,k,0) = fi(i,j,k) + 0.25 * ((1.0 - kk) * L1 + (1 + kk) * L2 * L) * A0;
                    FLR(i,j,k,1) = fi(i,j,k + 1) - 0.25 * ((1.0 - kk) * R1 + (1 + kk) * R2 * R) * A2;
                }
            }
}












// Array<double> WENO(int flag,Array<double> U,Array<double> V,Array<double> xnode,Array<double> ynode,Array<double> fi){
//     int ni = xnode.G
// }