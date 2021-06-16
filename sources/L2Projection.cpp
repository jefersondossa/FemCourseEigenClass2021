/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "L2Projection.h"
#include "PostProcess.h"
#include "tpanic.h"
#include <string.h>

L2Projection::L2Projection() {
}

L2Projection::L2Projection(int bctype, int materialid, MatrixDouble &proj, MatrixDouble Val1, MatrixDouble Val2) {
    projection = proj;
    BCType = bctype;
    BCVal1 = Val1;
    BCVal2 = Val2;
    this->SetMatID(materialid);
}

L2Projection::L2Projection(const L2Projection &copy) {
    projection = copy.projection;
    forceFunction = copy.forceFunction;
    SolutionExact = copy.SolutionExact;
    BCType = copy.BCType;
    BCVal1 = copy.BCVal1;
    BCVal2 = copy.BCVal2;

}

L2Projection &L2Projection::operator=(const L2Projection &copy) {
    projection = copy.projection;
    forceFunction = copy.forceFunction;
    SolutionExact = copy.SolutionExact;
    BCType = copy.BCType;
    BCVal1 = copy.BCVal1;
    BCVal2 = copy.BCVal2;
    return *this;
}

L2Projection *L2Projection::Clone() const {
    return new L2Projection(*this);
}

L2Projection::~L2Projection() {
}

MatrixDouble L2Projection::GetProjectionMatrix() const {
    return projection;
}

void L2Projection::SetProjectionMatrix(const MatrixDouble &proj) {
    projection = proj;
}

void L2Projection::Contribute(IntPointData &data, double weight, MatrixDouble &EK, MatrixDouble &EF) const {
    int nstate = this->NState();
    int nshape = data.phi.size();

    VecDouble result(data.x.size());
    MatrixDouble deriv(data.x.size(), data.x.size());

    // if (!SolutionExact)

    VecDouble val2(data.x.size());
    val2(0) =  BCVal2(0);
    if (SolutionExact) {
        SolutionExact(data.x, result, deriv);
        val2(0) = result(0);
    }
    VecDouble phi = data.phi;

    // std::cout << "BCVAL \n" << result<< std::endl;
    std::cout << "BCVAL2 \n" << val2<< std::endl;
    std::cout << "phi, weight \n" << phi << " " << weight << std::endl;

    

    switch (this->GetBCType()) {
        case 0:
        {
            EF += (MathStatement::gBigNumber * val2(0) * weight) * data.phi;
            EK += (MathStatement::gBigNumber * weight) * data.phi * data.phi.transpose();
            // for(auto iv = 0; iv < nstate; iv++){
			// 	for(auto in = 0 ; in < nshape; in++) {
			// 		EF(nstate*in+iv,0) += 100000000 * val2(iv) * phi(in) * weight;
			// 		for (auto jn = 0 ; jn < nshape; jn++) {
			// 			EK(nstate*in+iv,nstate*jn+iv) += phi(in) * phi(jn) * weight;
			// 		}//jn
			// 	}//in
			// }//iv
            // std::cout << "Matrix = \n" << EK << std::endl;
            // std::cout << "Rhs = \n" << EF << std::endl;
            break;
        }

        case 1:
        {
            for(auto iv = 0; iv < nstate; iv++){
				for(auto in = 0 ; in < nshape; in++) {
					EF(nstate*in+iv,0) += val2(iv) * phi(in) * weight;
				}//in
			}//iv
            break;
        }

        default:
        {
            std::cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " not implemented\n";
        }
    }
}

int L2Projection::NEvalErrors() const {
    return 3;
}

void L2Projection::ContributeError(IntPointData &data, VecDouble &u_exact, MatrixDouble &du_exact, VecDouble &errors) const {
    return;
}

int L2Projection::VariableIndex(const PostProcVar var) const {
    if (var == ESol) return ESol;
    if (var == EDSol) return EDSol;

    // Code should not reach this point. This return is only here to stop compiler warnings.
    DebugStop();
    return 0;
}

L2Projection::PostProcVar L2Projection::VariableIndex(const std::string & name) {
    if (!strcmp("Solution", name.c_str())) return ESol;
    if (!strcmp("Derivative", name.c_str())) return EDSol;

    // Code should not reach this point. This return is only here to stop compiler warnings.
    DebugStop();
    return ENone;
}

int L2Projection::NSolutionVariables(const PostProcVar var) {
    if (var == ESol) return this->NState();
    if (var == EDSol) return this->NState();

    // Code should not reach this point. This return is only here to stop compiler warnings.
    DebugStop();
    return 0;
}

void L2Projection::PostProcessSolution(const IntPointData &data, const int var, VecDouble &Solout) const {
    VecDouble sol = data.solution;
    int solsize = sol.size();
    int rows = data.dsoldx.rows();
    int cols = data.dsoldx.cols();
    MatrixDouble gradu(rows, cols);
    gradu = data.dsoldx;

    int nstate = this->NState();

    switch (var) {
        case 0: //None
        {
            std::cout << " Var index not implemented " << std::endl;
            // return;
            DebugStop();
        }
        case 1: //ESol
        {
            //+++++++++++++++++
            // Please implement me
            std::cout << "\nPLEASE IMPLEMENT ME\n" << __PRETTY_FUNCTION__ << std::endl;
            return;
            DebugStop();
            //+++++++++++++++++
        }
            break;

        case 2: //EDSol
        {
            //+++++++++++++++++
            // Please implement me
            std::cout << "\nPLEASE IMPLEMENT ME\n" << __PRETTY_FUNCTION__ << std::endl;
            return;
            DebugStop();
            //+++++++++++++++++
        }
            break;
        default:
        {
            std::cout << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
}
