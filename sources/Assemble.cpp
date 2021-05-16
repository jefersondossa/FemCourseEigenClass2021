/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Assemble.h"
#include "CompMesh.h"
#include "GeoMesh.h"
#include "MathStatement.h"
#include "CompElement.h"



Assemble::Assemble() {
}

Assemble::Assemble(CompMesh *mesh) {
    cmesh = mesh;
}

Assemble::Assemble(const Assemble &copy) {
    cmesh = copy.cmesh;
}

Assemble &Assemble::operator=(const Assemble &copy) {
    cmesh = copy.cmesh;
    return *this;
}

void Assemble::SetMesh(CompMesh *mesh) {
    cmesh = mesh;
}

int64_t Assemble::NEquations() {
    int64_t neq = 0;
    int64_t i, ncon = cmesh->GetDOFVec().size();
    for (i = 0; i < ncon; i++) {
        DOF dof = cmesh->GetDOF(i);
        int64_t dofsize = dof.GetNShape() * dof.GetNState();
        neq += dofsize;
    }
    return neq;
}

void Assemble::OptimizeBandwidth() {    
}

void Assemble::Compute(MatrixDouble &globmat, MatrixDouble &rhs) {
    int nelem = cmesh->GetGeoMesh()->NumElements();
    int ne = this->NEquations();

    int IG = 0, JG = 0;

    for (int el = 0; el < nelem; el++) {
        CompElement *cel = cmesh->GetElement(el);

        int nshape = cel->NShapeFunctions();
        int nstate = cel->GetStatement()->NState();
        MatrixDouble ek(nstate * nshape, nstate * nshape);
        MatrixDouble ef(nstate * nshape, 1);
        ek.setZero();
        ef.setZero();

        cel->CalcStiff(ek, ef);
//        
//        ek.Print();
//        ef.Print();

        // Implement the assembly method
        std::cout << "Please implement me\n";
        DebugStop();

    }
}
