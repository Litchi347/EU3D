







#include <fstream>
#include <iostream>
#include <mpi.h>





#include "Mesh.hpp"
#include "FileReader.hpp"
#include "Global.hpp"

using namespace std;














void Mesh::MeshProcess(char *gridname)
{

    FileReader grid;
    grid.readFile(gridname);
    total_ni = grid.getIntParameter("total_ni");
    xa = grid.getdoubleParameter("xa");
    xb = grid.getdoubleParameter("xb");

    total_nj = grid.getIntParameter("total_nj");
    ya = grid.getdoubleParameter("ya");
    yb = grid.getdoubleParameter("yb");

    total_nk = grid.getIntParameter("total_nk");
    za = grid.getdoubleParameter("za");
    zb = grid.getdoubleParameter("zb");

    bc = grid.getIntParameter("bc");


    Initial();


    Boundary();



}











void Mesh::Initial()
{
    dx = (xb - xa) / (total_ni - 1);
    ni = (total_ni - 1) / m_block_x + 1;
    xnode.Initial(ni + 2 * bc);
    for (int i = bc; i < ni + bc; i++)
        xnode(i) = xa + (myid_x * (ni - 1) + i - bc) * dx;

    dy = (yb - ya) / (total_nk - 1);
    nj = (total_nj - 1) / m_block_y + 1;
    ynode.Initial(nj + 2 * bc);
    for (int j = bc; j < nj + bc; j++)
        ynode(j) = ya + (myid_y * (nk - 1) + j - bc) * dy;

    dz = (zb - za) / (total_nk - 1);
    nk = (total_nk - 1) / m_block_z + 1;
    znode.Initial(nk + 2 * bc);
    for (int k = bc; k < nk + bc; k++)
        znode(k) = za + (myid_z * (nk - 1) + k - bc) * dz;

    total_xnode.Initial(total_ni + 2 * bc);
    total_ynode.Initial(total_nj + 2 * bc);
    total_znode.Initial(total_nk + 2 * bc);
}











void Mesh::Boundary()
{
    for (int i = 0; i < bc; i++)
    {
        xnode(i) = 2 * xnode(bc) - xnode(2 * bc - i);
        xnode(ni + bc + i) = 2 * xnode(ni + bc - 1) - xnode(ni + bc - 1 - i - 1);
        ynode(i) = 2 * ynode(bc) - ynode(2 * bc - i);
        ynode(nj + bc + i) = 2 * ynode(nj + bc - 1) - ynode(nj + bc - 1 - i - 1);
        znode(i) = 2 * znode(bc) - znode(2 * bc - i);
        znode(nk + bc + i) = 2 * znode(nk + bc - 1) - znode(nk + bc - 1 - i - 1);
    }
}











void Mesh::MeshCheck()
{
    string filename = "./output_" + to_string(numprocs) + "/grid/grid_check_" + to_string(myid) + ".dat";
    std::ofstream outfile(filename);
    outfile << "VARIABLEs = X,Y,Z" << '\n';
    outfile << "ZONE I=" << ni + 2 * bc << '\t' << "J=" << nj + 2 * bc << '\t' << "K=" << nk + 2 * bc << '\n';
    outfile << "datapacking=block" << '\n';

    for (int k = 0; k < nk + 2 * bc; k++)
        for(int j = 0; j < nj + 2 * bc; j++)
            for(int i = 0; i < ni + 2 * bc; i++)
                outfile << xnode(i) << '\n';
    for (int k = 0; k < nk + 2 * bc; k++)
		for (int j = 0; j < nj + 2 * bc; j++)
			for (int i = 0; i < ni + 2 * bc; i++)
				outfile << ynode(j) << '\n';
    for (int k = 0; k < nk + 2 * bc; k++)
		for (int j = 0; j < nj + 2 * bc; j++)
			for (int i = 0; i < ni + 2 * bc; i++)
				outfile << znode(k) << '\n';
    outfile.close();
}











Array<double, 1> Mesh::GetXnode()
{
    return xnode;
}

Array<double, 1> Mesh::GetYnode()
{
    return ynode;
}