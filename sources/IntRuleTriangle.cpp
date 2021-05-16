/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <iostream> 
#include "IntRuleTriangle.h"
#include "tpanic.h"

IntRuleTriangle::IntRuleTriangle(){
    
}

IntRuleTriangle::IntRuleTriangle(int order) {
    SetOrder(order);
}

void IntRuleTriangle::SetOrder(int order) {
    if ((order < 0) || (order > MaxOrder())) DebugStop();
    fOrder = order;
    
    switch (order) {
        case 0:
        case 1:
            fPoints.resize(1,Dimension());
            fWeights.resize(1);
            fPoints(0,0) = 1./3.;
            fPoints(0,1) = 1./3.;
            fWeights(0) = 1.;
        case 2:
            fPoints.resize(3,Dimension());
            fWeights.resize(3);
            fPoints(0,0) = 0.5;
            fPoints(0,1) = 0.;
            fPoints(1,0) = 0.5;
            fPoints(1,1) = 0.5;
            fPoints(2,0) = 0.;
            fPoints(2,1) = 0.5;
            fWeights(0) = 1./3.;
            fWeights(1) = 1./3.;
            fWeights(2) = 1./3.;
        case 3:
            fPoints.resize(4,Dimension());
            fWeights.resize(4);
            fPoints(0,0) = 1./3.;
            fPoints(0,1) = 1./3.;
            fPoints(1,0) = 2./15.;
            fPoints(1,1) = 2./15.;
            fPoints(2,0) = 11./15.;
            fPoints(2,1) = 2./15.;
            fPoints(3,0) = 2./15.;
            fPoints(3,1) = 11./15.;
            fWeights(0) = -27./48.;
            fWeights(1) = 25./48.;
            fWeights(2) = 25./48.;
            fWeights(3) = 25./48.;
        case 4:
        case 5:
            fPoints.resize(7,Dimension());
            fWeights.resize(7);
            fPoints(0,0) = 1./3.;
            fPoints(0,1) = 1./3.;
            fPoints(1,0) = 1. - 2. * (6. - sqrt(15.)) / 21.;
            fPoints(1,1) = (6. - sqrt(15.)) / 21.;
            fPoints(2,0) = (6. - sqrt(15.)) / 21.;
            fPoints(2,1) = 1. - 2. * (6. - sqrt(15.)) / 21.;
            fPoints(3,0) = (6. - sqrt(15.)) / 21.;
            fPoints(3,1) = (6. - sqrt(15.)) / 21.;
            fPoints(4,0) = (6. + sqrt(15.)) / 21.;
            fPoints(4,1) = (6. + sqrt(15.)) / 21.;
            fPoints(5,0) = 1. - 2. * (6. + sqrt(15.)) / 21.;
            fPoints(5,1) = (6. + sqrt(15.)) / 21.;
            fPoints(6,0) = (6. + sqrt(15.)) / 21.;
            fPoints(6,1) = 1. - 2. * (6. + sqrt(15.)) / 21.;
            fWeights(0) = 0.1125;
            fWeights(1) = (155. - sqrt(15.)) / 2400.;
            fWeights(2) = (155. - sqrt(15.)) / 2400.;
            fWeights(3) = (155. - sqrt(15.)) / 2400.;
            fWeights(4) = (155. + sqrt(15.)) / 2400.;
            fWeights(5) = (155. + sqrt(15.)) / 2400.;
            fWeights(6) = (155. + sqrt(15.)) / 2400.;
        break;
    default:
        break;
    }
}
