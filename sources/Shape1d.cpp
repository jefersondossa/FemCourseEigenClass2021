//
//  Shape1d.h
//  FemSC
//
//  Created by Philippe Devloo on 03/04/18.
//

#include <cmath>
#include <math.h>
#include "tpanic.h"
#include "DataTypes.h"
#include "Shape1d.h"

void Shape1d::Shape(const VecDouble &xi, VecInt &orders, VecDouble &phi, MatrixDouble &dphi){
    
    if (orders[0] < 0 || orders[1] < 0 || orders[2] < 0) {
        std::cout << "Shape1d::Shape: Invalid dimension for arguments: order\n";
        DebugStop();
    }
    if (orders[0] > 1 || orders[1] > 1) {
        std::cout << "Shape1d::Shape: Invalid dimension for arguments: order\n";
        DebugStop();
    }
    if (orders[2] > 2) {
        std::cout << "Shape1d::Shape: Please implement it for order > 2\n";
        DebugStop();
    }
    
    auto nshape = NShapeFunctions(orders);
    phi.resize(nshape);
    dphi.resize(1,nshape);
    phi.setZero(); dphi.setZero();

    phi[0] = 0.5 * (1. - xi(0));
    phi[1] = 0.5 * (1. + xi(0));
    dphi(0,0) = -0.5;
    dphi(0,1) =  0.5;

    int count = 2;
    int is;
    for (is = 2; is < 3; is++) {
        if(orders[is] == 2){
            int is1 = SideNodeLocIndex(is, 0);
            int is2 = SideNodeLocIndex(is, 1);
            phi[is] = 4.*phi[is1] * phi[is2];
            dphi(0, is) = 4.*(dphi(0, is1) * phi[is2] + phi[is1] * dphi(0, is2));
            count++;
        } else if (orders[is] != 1) DebugStop();
    }

    if(count != nshape) DebugStop();
    for(int is = 3 ; is< nSides; is++) if(orders[is] != 1 && orders[is] != 2) DebugStop();
}

/// returns the number of shape functions associated with a side
int Shape1d::NShapeFunctions(int side, int order){

    if(order < 1 || order >2) DebugStop();
    switch (side)
    {
    case 0:
        return 1;
        break;
    case 1:
        return 1;
        break;
    case 2:
        return order-1;
        break;
    
    default:
        std::cout << "Shape1d::NShapeFunctions : Wrong side " << side << "\n";
        DebugStop();
        return -1;
        break;
    }
    return -1;
}

/// returns the total number of shape functions
int Shape1d::NShapeFunctions(VecInt &orders) {
    
    int nsf_tot = 0;
    for (int is=0; is<3; is++) {
        nsf_tot += NShapeFunctions(is, orders[is]);
    }
    
    return nsf_tot;
}
