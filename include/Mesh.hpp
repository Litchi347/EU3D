// 在并行环境下的三维结构化网格生成与管理




#ifndef MESH_H
#define MESH_H





# include <algorithm>
# include <iostream>




# include "Array.hpp"

using namespace std;
using namespace ARRAY;

class Mesh
{
private:

    int total_ni;
    Array<double, 1> total_xnode;
    int ni;
    double xa;
    double xb;
    Array<double, 1> xnode;
    double dx;


    int total_nj;
    Array<double, 1> total_ynode;
    int nj;
    double ya;
    double yb;
    Array<double, 1> ynode;
    double dy;


    int total_nk;
    Array<double, 1> total_znode;
    int nk;
    double za;
    double zb;
    Array<double, 1> znode;
    double dz;

    int bc;
    
public:
    friend class Euler;


    Mesh() = default;


    void MeshProcess(char *gridname);


    void Initial();


    void Boundary();


    void MeshCheck();


    ~Mesh(){ ; };


    Array<double,1> GetXnode();
    Array<double,1> GetYnode();
};
#endif