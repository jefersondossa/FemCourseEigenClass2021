/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "GeomQuad.h"

GeomQuad::GeomQuad() {
}

GeomQuad::~GeomQuad() {
}

GeomQuad::GeomQuad(const GeomQuad &copy) {
    fNodeIndices = copy.fNodeIndices;
}

GeomQuad& GeomQuad::operator=(const GeomQuad& copy) {
    fNodeIndices = copy.fNodeIndices;
    return *this;
}

void GeomQuad::Shape(const VecDouble &xi, VecDouble &phi, MatrixDouble &dphi) {
    std::cout << "\nPLEASE IMPLEMENT ME\n" << __PRETTY_FUNCTION__ << std::endl;
    if(xi.size() != Dimension || phi.size() != nCorners || dphi.rows() != Dimension || dphi.cols() != nCorners) DebugStop();

    double qsi = xi[0];
    double eta = xi[1];

    phi[0] = 0.25 * (1.-qsi) * (1.-eta);
    phi[1] = 0.25 * (1.+qsi) * (1.-eta);
    phi[2] = 0.25 * (1.+qsi) * (1.+eta);
    phi[3] = 0.25 * (1.-qsi) * (1.+eta);

    dphi(0, 0) = -0.25 * (1.-eta);
    dphi(1, 0) = -0.25 * (1.-qsi);

    dphi(0, 1) =  0.25 * (1.-eta);
    dphi(1, 1) = -0.25 * (1.+qsi);

    dphi(0, 2) =  0.25 * (1.+eta);
    dphi(1, 2) =  0.25 * (1.+qsi);

    dphi(0, 3) = -0.25 * (1.+eta);
    dphi(1, 3) =  0.25 * (1.-qsi);

}

void GeomQuad::X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x) {
    std::cout << "\nPLEASE IMPLEMENT ME\n" << __PRETTY_FUNCTION__ << std::endl;
    if(xi.size() != Dimension) DebugStop();
    if(x.size() < NodeCo.rows()) DebugStop();
    if(NodeCo.cols() != nCorners) DebugStop();
    VecDouble phi(nCorners);
    MatrixDouble dphi(Dimension, nCorners);
    
    Shape(xi, phi, dphi);
    int space = NodeCo.rows();

    for (int i = 0; i < Dimension; i++) {
        for (int j = 0; j < nCorners; j++) {
            x[i] += phi[j] * NodeCo(i, j);
        }
    }   
}

void GeomQuad::GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, MatrixDouble &gradx) {
    if(xi.size() != Dimension) DebugStop();
    if(x.size() < NodeCo.rows()) DebugStop();
    if(NodeCo.cols() != nCorners) DebugStop();
    gradx.setZero();
    x.setZero();
    int nrow = NodeCo.rows();
    int ncol = NodeCo.cols();

    VecDouble phi(nCorners);
    MatrixDouble dphi(Dimension, nCorners);
    Shape(xi, phi, dphi);
    
    for (int i = 0; i < nCorners; i++) {
        for (int j = 0; j < Dimension; j++) {
            x[j] += NodeCo(j,i) * phi[i];
            gradx(j,0) += NodeCo(j,i) * dphi(0,i);
            gradx(j,1) += NodeCo(j,i) * dphi(1,i);
        }
    }
}

void GeomQuad::SetNodes(const VecInt &nodes) {
    if(nodes.size() != nCorners) DebugStop();
    fNodeIndices = nodes;
}

void GeomQuad::GetNodes(VecInt &nodes) const{
    nodes = fNodeIndices;
}

int GeomQuad::NodeIndex(int node) const {
    return fNodeIndices[node];
}

int GeomQuad::NumNodes() {
    return nCorners;
}

GeoElementSide GeomQuad::Neighbour(int side) const {
    return fNeighbours[side];
}

void GeomQuad::SetNeighbour(int side, const GeoElementSide &neighbour) {
    fNeighbours[side] = neighbour;
}
