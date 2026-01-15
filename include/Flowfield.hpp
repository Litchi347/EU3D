





#ifndef FlowField
#define FlowField





# include <iostream>





# include"Array.hpp"
# include "TimeAdv.hpp"
# include "Reaction.hpp"
# include "Function.hpp"

using namespace std;
using namespace ARRAY;

class Flowfield
{
private:

    int ni;
    Array<double,1> xnode;
    double dx;

    int nj;
    Array<double,1> ynode;
    double dy;

    int nk;
    Array<double,1> znode;
    double dz;

    int bc;

    double dt;
    double time;

    double ft1 = 0, ft2 = 0, ft3 = 0, ft4 = 0, ft5 = 0, ft6 = 0;

    int NS;
    int NR;
    Reaction React;
    const double R = 8.31434;
    Array<double,1> Mw;
    Array<double,1> Ri;
    Array<double,2> Coeff0;
    Array<double,2> Coeff1;

    TimeAdv Time;


    double M_wave;
    double U_post;
    double P_pre;
    double P_post;
    double T_pre;
    double T_post;

    double cita;
    double V3;
    Array<double,3> U;
    Array<double,3> V;
    Array<double,3> W;
    Array<double,3> P;
    double P_bound;
    Array<double,3> D;
    Array<double,3> T;
    double T_bound;
    Array<double,3> C;
    Array<double,3> Ma;
    Array<double,3> Wav;
    Array<double,3> Rgas;
    Array<double,3> Cp;
    Array<double,3> H;
    Array<double,3> E;
    Array<double,3> Gamma;
    Array<double,4> Mr;
    Array<double,4> Mc;
    Array<double,4> Mi;
    Array<double,4> Yi;
    Array<double,4> Di;
    Array<double,1> Cpi;
    Array<double,1> Hi;
    Array<double,1> Ei;

    Array<double,4> F;
    Array<double,4> G;
    Array<double,4> Q;
    Array<double,4> CS;

    Function Fun;


    Array<double,2> Uint, Vint, Wint, Pint, Dint, Tint, Hint, Eint, Gint;
    Array<double,3> Yint;


    Array<double,1> Yi_temp0, YL, YR;
    Array<double,4> PLR, DLR, ULR, VLR, WLR, HLR, GLR, YLR, YL_temp, YR_temp;
    Array<double,3> Yi_temp;
    Array<double,4> Partion_T;
    Array<double,4> RHS;


    Array<double,1> m_U_s, m_V_s, m_P_s, m_D_s, m_C_s, m_Gamma_s, m_H_s, m_T_s,m_Yi_s;
    Array<double,1> m_U_r, m_V_r, m_P_r, m_D_r, m_C_r, m_Gamma_r, m_H_r, m_T_r,m_Yi_r;
    Array<double,1> send_data_1, recv_data_1, send_data_2, recv_data_2, send_data_3,recv_data_3;

    struct SendRecv_Data
    {
        Array<double,1> D, U, V, P, T, H, Gamma, C, Yi;
        void Initial(int size,int _NS)
        {
            D.Initial(size);
            U.Initial(size);
            V.Initial(size);
            P.Initial(size);
            T.Initial(size);
            H.Initial(size);
            Gamma.Initial(size);
            C.Initial(size);
            Yi.Initial(size * _NS);
        }
    };
    SendRecv_Data send_data, recv_data;

public:
    friend class Euler;


    Flowfield() = default;


    void InputRead(char *initialization, Array<double,1> &xnode,Array<double,1> &ynode,Array<double,1> &znode, int bc);


    void Construction();


    void ReConstruction(int meshnum);


    void FieldInitial(Array<double,1> &Ri, Array<double,1> &Mw, Array<double,4> &Mi_temp, Array<double,4> &Yi_temp,Array<double,2> &coeff0, Array<double,2> &coeff1);
    
    
    void CFLcondition(double cfl, double Final_Time);


    void FieldBoundary_3dODM();
    void FieldBoundary_1dZND();
    void FieldBoundary_3dRSBI();


    void Advection(int, int);


    void Update_after_Adv();


    void Explicit();
    void Explicit(int step, int i, int j, int k);


    void Update_IMEX(Array<double, 4> &Wi, Array<double, 5> &MD);


    void Update_after_CS();
    void Update_after_CS(int step, int i, int j, int k);


    ~Flowfield(){ ; };


    void AUSM(int, Array<double,4> &, void (*Diff)(int, Array<double,1> &, Array<double,1> &, Array<double,1> &, Array<double,4> &, Array<double,3> &, int ,int));


    double Get_temp(double, int,int ,int);
    double Get_temp(double, int,int ,int,Array<double,1> &);


    void GetPartial_T();


    void Mpi_Boundary();




    void PackagePrev(std::vector<double> &, int i, int j, int k);
    void PackageUpdate(std::vector<double> &,int meshnum);
    void PackageUpdate(std::vector<double> &,int, std::vector<int> &, std::vector<int> &);
    void PackageUpdate(std::vector<double> &, std::vector<int> &, std::vector<int> &);


    void UnpackagePrev(std::vector<double> &, int meshnum);
    void UnpackagePrev(std::vector<double> &, int, std::vector<int> &, std::vector<int> &);
    void UnpackagePrev(std::vector<double> &, std::vector<int> &, std::vector<int> &);
    void UnpackagePrev(std::vector<double> &, Array<int, 1> &, Array<int,1> &);
    void UnpackageUpdate(std::vector<double> &, std::vector<int> &, int, std::vector<int> &, std::vector<int> &);
    void UnpackageUpdate(std::vector<double> &, std::vector<int> &, std::vector<int> &, std::vector<int> &);
};
#endif