

#include <iostream>
#include <math.h>
#include <fstream>

#include "IntRule.h"
#include "IntRule1d.h"
#include "IntRuleQuad.h"
#include "IntRuleTetrahedron.h"
#include "IntRuleTriangle.h"
#include "Topology1d.h"
#include "TopologyTriangle.h"
#include "TopologyQuad.h"
#include "TopologyTetrahedron.h"
#include "DataTypes.h"
#include "Analysis.h"
#include "VTKGeoMesh.h"
#include "Poisson.h"
#include "L2Projection.h"
#include "PostProcess.h"
#include "CompMesh.h"
#include "ReadGmsh.h"

using std::cout;
using std::endl;
using std::cin;

auto exactSol = [](const VecDouble &loc,
  VecDouble &u,
  MatrixDouble &gradU){
  const auto &x=loc[0];
  const auto &y=loc[1];
  
  // const auto &d = 1.; // distance betweel injection and production wells
  u[0]= x*x-y*y;//log(hypot(x,y)) - log(hypot(x-d,y-d)) - log(hypot(x+d,y-d)) - log(hypot(x-d,y+d)) - log(hypot(x+d,y+d));
  gradU(0,0) = 2.*x;//x/(x*x+y*y) - (x-d)/(pow(x-d,2)+pow(y-d,2)) - (x+d)/(pow(x+d,2)+pow(y-d,2)) - (x-d)/(pow(x-d,2)+pow(y+d,2)) - (x+d)/(pow(x+d,2)+pow(y+d,2));
  gradU(1,0) = -2.*y;//y/(x*x+y*y) - (y-d)/(pow(x-d,2)+pow(y-d,2)) - (y-d)/(pow(x+d,2)+pow(y-d,2)) - (y+d)/(pow(x-d,2)+pow(y+d,2)) - (y+d)/(pow(x+d,2)+pow(y+d,2));
  // gradU(2,0) = 0;//optional
};

int main (){

    ReadGmsh *reader;
    reader = new ReadGmsh();
    GeoMesh gmesh;
    reader -> Read(gmesh,"../meshes/1element.msh");
    std::stringstream text_name;
    text_name << "geometry.txt";
    std::ofstream textfile(text_name.str().c_str());
    gmesh.Print(textfile);
    VTKGeoMesh::PrintGMeshVTK(&gmesh, "geometry.vtk");

    CompMesh cmesh(&gmesh);
    cmesh.SetDefaultOrder(1);
    MatrixDouble perm(3,3);
    perm.setZero();
    perm(0,0) = 1.;
    perm(1,1) = 1.;
    perm(2,2) = 1.;
    Poisson *mat1 = new Poisson(1,perm);
    MatrixDouble proj(1,1), val1(1,1), val2(1,1), val3(1,1);
    proj.setZero();
    val1.setZero();
    val2.setZero();
    val1(0,0) = 1.;
    val3(0,0) = 1.;
    L2Projection *bc_bottom = new L2Projection(0,2,proj,val1,val2);
    L2Projection *bc_right = new L2Projection(0,3,proj,val1,val2);
    L2Projection *bc_top = new L2Projection(0,4,proj,val1,val2);
    L2Projection *bc_left = new L2Projection(0,5,proj,val1,val3);
    // L2Projection *bc_point = new L2Projection(0,6,proj,val1,val2);
    bc_bottom->SetExactSolution(exactSol);
    bc_right->SetExactSolution(exactSol);
    bc_top->SetExactSolution(exactSol);
    bc_left->SetExactSolution(exactSol);
    // bc_point->SetExactSolution(exactSol);
    std::vector<MathStatement *> mathvec = {0,mat1,bc_bottom,bc_right,bc_top,bc_left};


    cmesh.SetMathVec(mathvec);

    cmesh.AutoBuild();
    cmesh.Resequence();

    std::stringstream cmesh_name;
    cmesh_name << "cmesh.txt";
    std::ofstream textfile2(cmesh_name.str().c_str());
    cmesh.Print(cmesh_name);
    VTKGeoMesh::PrintCMeshVTK(&cmesh,2, "cmesh.vtk");


    Analysis Analysis(&cmesh);
    Analysis.RunSimulation();

    // std::cout << "Matrix = \n" << Analysis.GlobalSystem << std::endl;
    // std::cout << "Rhs = \n" << Analysis.RightHandSide << std::endl;
    // std::cout << "Solution = \n" << Analysis.Solution << std::endl;

    VTKGeoMesh::PrintCMeshVTK(&cmesh,2, "cmesh.vtk");
    PostProcess *post;
    std::vector<std::string> scalarVars(1), vectorVars(1);
    scalarVars[0] = "ESol";
    vectorVars[0] = "EDSol";
    // Analysis.DefineGraphMesh(dim,scalarVars,vectorVars,"SolutionH1.vtk");
    // constexpr int resolution{0};
    // Analysis.PostProcess(resolution,dim);	
    // post = new PostProcess(&Analysis);
    // Analysis.PostProcessSolution("sol.txt",post);
    return 0;
}
