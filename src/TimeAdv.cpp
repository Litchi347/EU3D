







#include <iostream>





#include "TimeAdv.hpp"













void TimeAdv::EE(int NS, double dt, Array<double, 1> &xnode, Array<double, 1> &ynode, Array<double, 1> &znode, int bc, Array<double, 4> &F, Array<double, 4> &G, Array<double, 4> &Q, Array<double, 4> &RHS)
{

    int ni = xnode.GetSize() - 2 * bc;
    int nj = ynode.GetSize() - 2 * bc;
    int nk = znode.GetSize() - 2 * bc;

    for (int i = bc; i < ni + bc; i++)
        for (int j = bc; j < nj + bc; j++)
            for (int k = bc; k < nk + bc; k++)
                for (int s = 0; s < NS + 4; s++)
                    RHS(i, j, k, s) = -dt / (xnode(i) - xnode(i - 1)) * (F(i, j, k, s) - F(i - 1, j, k, s)) -
                                      dt / (ynode(j) - ynode(j - 1)) * (G(i, j, k, s) - G(i, j - 1, k, s)) -
                                      dt / (znode(k) - znode(k - 1)) * (Q(i, j, k, s) - G(i, j, k - 1, s));
}













void TimeAdv::TVD_RK3(int NS, double dt, Array<double, 1> &xnode, Array<double, 1> &ynode, Array<double, 1> &znode, int bc, Array<double, 4> &F, Array<double, 4> &G, Array<double, 4> &Q, Array<double, 4> &RHS)
{
    int ni = xnode.GetSize();
    int nj = ynode.GetSize();
    int nk = znode.GetSize() - 2 * bc;
    for (int i = bc; i < ni + bc; i++)
        for (int j = bc; j < nj + bc; j++)
            for (int k = bc; k < nk + bc; k++)
                for (int s = 0; s < NS + 4; s++)
                    RHS(i, j, k, s) = RHS(i, j, k, s);
}