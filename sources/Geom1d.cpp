/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "Geom1d.h"

Geom1d::Geom1d() {
}

Geom1d::~Geom1d() {
}

Geom1d::Geom1d(const Geom1d &copy) {
    fNodeIndices = copy.fNodeIndices;
}

Geom1d& Geom1d::operator=(const Geom1d& copy) {
    fNodeIndices = copy.fNodeIndices;
    return *this;
}

void Geom1d::Shape(const VecDouble &xi, VecDouble &phi, MatrixDouble &dphi) {
    if (xi.size() == 0) DebugStop();
    phi.resize(2);
    dphi.resize(1,2);
    phi(0) = 0.5 * (1. - xi(0));
    phi(1) = 0.5 * (1. + xi(0));
    dphi(0,0) = -0.5;
    dphi(0,1) =  0.5;
}

void Geom1d::X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x) {
    if (x.size() == 0) DebugStop();
    VecDouble phi;
    MatrixDouble dphi;
    Shape(xi,phi,dphi);
    for (int i = 0; i < NodeCo.cols(); i++){
        x(0) += NodeCo(0,i) * phi(i);
    }
    
}

void Geom1d::GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, MatrixDouble &gradx) {
    if (x.size() == 0) DebugStop();
    VecDouble phi;
    MatrixDouble dphi;
    Shape(xi,phi,dphi);
    for (int i = 0; i < NodeCo.cols(); i++){
        gradx(0,0) += NodeCo(0,i) * dphi(0,i);
    }
}

void Geom1d::SetNodes(const VecInt &nodes) {
    if(nodes.rows() != 2)
    fNodeIndices = nodes;
}

void Geom1d::GetNodes(VecInt &nodes) const{
    nodes = fNodeIndices;
}

int Geom1d::NodeIndex(int node) const{
    return fNodeIndices[node];
}

int Geom1d::NumNodes() {
    return nCorners;    
}

GeoElementSide Geom1d::Neighbour(int side) const {
    return fNeighbours[side];
}

void Geom1d::SetNeighbour(int side, const GeoElementSide &neighbour) {
    fNeighbours[side]=neighbour;
}
