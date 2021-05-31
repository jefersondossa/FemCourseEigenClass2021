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
    phi(0) = 0.5 * (1. - xi(0));
    phi(1) = 0.5 * (1. + xi(0));
    dphi(0,0) = -0.5;
    dphi(0,1) =  0.5;
}

void Geom1d::X(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x) {
    if(xi.size() != Dimension) DebugStop();
    // if(x.size() != NodeCo.rows()) DebugStop();
    if(NodeCo.cols() != nCorners) DebugStop();
    VecDouble phi(nCorners);
    MatrixDouble dphi(Dimension, nCorners);

    x.setZero();
    Shape(xi, phi, dphi);
    int nrow = NodeCo.rows();
    int ncol = NodeCo.cols();

    for (int i = 0; i < nrow; i++) {
        x[i] = 0.0;
        for (int j = 0; j < ncol; j++) {
            x[i] += phi[j] * NodeCo(i, j);
        }
    }
}

void Geom1d::GradX(const VecDouble &xi, MatrixDouble &NodeCo, VecDouble &x, MatrixDouble &gradx) {
    if(xi.size() != Dimension) DebugStop();
    if(x.size() != NodeCo.rows()) DebugStop();
    if(NodeCo.cols() != nCorners) DebugStop();
    gradx.resize(Dimension, 1);
    gradx.setZero();
    x.resize(Dimension);
    x.setZero();
    int nrow = NodeCo.rows();
    int ncol = NodeCo.cols();

    VecDouble phi(nCorners);
    MatrixDouble dphi(Dimension, nCorners);
    Shape(xi, phi, dphi);
    for (int i = 0; i < ncol; i++) {
        for (int j = 0; j < nrow; j++) {
            x[j] += NodeCo(j,i) * phi[i];
            gradx(j, 0) += NodeCo(j, i) * dphi(0, i);
        }
    }
}

void Geom1d::SetNodes(const VecInt &nodes) {
    if(nodes.rows() != 2) DebugStop();
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
