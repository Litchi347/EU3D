




#ifndef EULER_H
#define EULER_H





# include <fstream>
# include <iostream>
# include <chrono>





# include "Array.hpp"
# include "Mesh.hpp"
# include "Flowfield.hpp"
# include "Reaction.hpp"
# include "Global.hpp"
# include "LoadBalance.hpp"

using namespace ARRAY;
class Euler
{
private:
    char ctrlfile[100] = "./input/ctrl.txt";
    char meshFile[100];
    char inputFile[100];
    char Reaction_Model[100];
    char thermoFile[100];
    char react[50] = {0}, diff[50] = {0}, timeadv[50] = {0};
    int Reaction_Sch;
    int type;
    int Diff_Sch;
    int TimeAdv_Sch;
    double Final_Time;
    double cfl;
    double dt;
    int count;
    int iteration;

    double tr1 = 0, tr2 = 0, tr3 = 0, tr4 = 0;
    double dt1 = 0.0, dt2 = 0.0, dt3 = 0.0, dt4 = 0.0, dt5 = 0.0, dt6 = 0.0, dt7 = 0.0, dt8 = 0.0;
    double trackt = 0.0;

    Mesh Mymesh;
    Flowfield TwoDim;
    Reaction React;


    Flowfield ExtraFlow;
    Reaction ExtraReact;
    DynamicLoadBalancer DLB;
    int DLB_step = 0;
    double DLB_tol = 0;
    int pos = 0;
    int TransferNchem = 0;
    int AverageNchem = 0;
    Array<int,1> NchemIndex;
    std::vector<int> TransferIndex;
    int SendTransferMeshNum = 0;
    int RecvTransferMeshNum = 0;
    std::vector<double> PrevFlowData;
    std::vector<double> PrevFlowDataLocal;
    std::vector<double> SendPrevFlowDataLocal;
    std::vector<double> RecvPrevFlowDataLocal;
    std::vector<double> UpdateFlowData;
    std::vector<double> UpdateFlowDataLocal;
    std::vector<double> SendUpdateFlowDataLocal;
    std::vector<double> RecvUpdateFlowDataLocal;
    std::vector<double> PrevFlowNumID;
    std::vector<double> UpdateFlowNumID;
    Array<int,1> SendRecvNchem;
    Array<int,1> SendRecvMeshNumLocal;
    Array<int,1> SendRecvMeshNum;
    Array<int,1> SendRecvMeshNumTemp;
    Array<int,1> RecvPrevFlowMeshNum;
    Array<int,1> GatherRecvCounts;
    Array<int,1> GatherDispls;
    int sendid = 0;
    int recvid = 0;
    Array<int,1> Capacity;
    int Capacity_local = 0;
    MPI_Request request1[2];
    MPI_Status statuses1[2];
    MPI_Request request2[2];
    MPI_Status statuses2[2];
    
    
    void
    FileRead();
    
public:

    Euler();


    void Mpi_Initial();


    void Computing(int count);


    void Trapezoid(int type);


    void IMEX();


    void DNN(int type);


    void Output(int num);


    void Output_Total_Time(double t1, double t2, double t3, double t4, double t5, double t6, double t7, double t8, double t9, double t10);
    void Output_Compute_Time(double t1, double t2, double t3, double t4, double t5);

    ~Euler()
    {
        ;
    };
};
#endif