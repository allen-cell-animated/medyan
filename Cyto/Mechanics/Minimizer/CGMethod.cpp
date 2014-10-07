//
//  CGMethod.cpp
//  Cyto
//
//  Created by James Komianos on 9/10/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "CGMethod.h"
#include "ForceFieldManager.h"
#include <cmath>

using namespace std;

const double phi = (1 + sqrt(5)) / 2;
const double resphi = 2 - phi;

double CGMethod::GradSquare()
{
    double g = 0;
	for(auto it: *BeadDB::Instance(getBeadDBKey())) {
        
        g += (*it).CalcForceSquare();
	}
    
    return g;
}

double CGMethod::GradSquare(int i)
{
    double g = 0;
	for(auto it: *BeadDB::Instance(getBeadDBKey())) {
        
        g += (*it).CalcForceSquare(i);
	}
    
    return g;
}

double CGMethod::GradDotProduct()
{
    double g = 0;
	for(auto it: *BeadDB::Instance(getBeadDBKey())) {
        
        g += (*it).CalcDotForceProduct();
	}
    
    return g;
}


void CGMethod::MoveBeads(double d)
{
	for(auto it: *BeadDB::Instance(getBeadDBKey())) {
        
        (*it).coordinate[0] = (*it).coordinate[0] + d* (*it).force[0];
        (*it).coordinate[1] = (*it).coordinate[1] + d* (*it).force[1];
        (*it).coordinate[2] = (*it).coordinate[2] + d* (*it).force[2];
        
	}
    
    
}

void CGMethod::ShiftGradient(double d)
{
	for(auto it: *BeadDB::Instance(getBeadDBKey())) {
        
        (*it).force[0] = (*it).forceAux[0] + d* (*it).force[0];
        (*it).force[1] = (*it).forceAux[1] + d* (*it).force[1];
        (*it).force[2] = (*it).forceAux[2] + d* (*it).force[2];
	}
}




void CGMethod::PrintForces()
{
	cout << "Print Forces" << endl;
    for(auto it: *BeadDB::Instance(getBeadDBKey())) {
        
		for (int i = 0; i<3; i++) { cout << (*it).coordinate[i] << "  "<< (*it).force[i]<<"  "<<(*it).forceAux[i]<<endl;}
	}

    cout << "End of Print Forces" << endl;
}



double CGMethod::GoldenSection1(ForceFieldManager& FFM)
{
	
    const double EPS = 1e-6;
	double a = 0;
	double b = 100;
    double phi = 0.5 * (1 + sqrt(5) );
    double inv_phi = 1/phi;
	double x1 = b - inv_phi * (b - a);
	double x2 = a + inv_phi * (b - a);
    
    
    
	while (abs(b - a) > EPS)
	{
		if (FFM.ComputeEnergy(x1) >= FFM.ComputeEnergy(x2) ){
            a = x1;
            x1 = x2;
            x2 =a + inv_phi * (b - a);
        }
        else {
            b = x2;
            x2 = x1;
            x1 = b - inv_phi * (b - a);
        }
    }
    
    
	return (a + b)/2.0;
}

double CGMethod::GoldenSection2(ForceFieldManager& FFM, double a, double b, double c, double tau)
{
    // a and c are the current bounds; the minimum is between them.
    // b is a center point
    // f(x) is some mathematical function elsewhere defined
    // a corresponds to x1; b corresponds to x2; c corresponds to x3
    // x corresponds to x4
    // tau is a tolerance parameter; see above
    double x;
    if (c - b > b - a)
        x = b + resphi * (c - b);
    else
        x = b - resphi * (b - a);
    if (abs(c - a) < tau * (abs(b) + abs(x)))
        return (c + a) / 2;
    
    double fx = FFM.ComputeEnergy(x);
    double fb = FFM.ComputeEnergy(b);
    
//    std::cout << "x = " << x << " b = " << b << std::endl;
//    std::cout << "fx = " << fx << " fb = " << fb << std::endl;
    
   // assert(fx != fb);
    
    if (fx < fb) {
        if (c - b > b - a) return GoldenSection2(FFM, b, x, c, tau);
        else return GoldenSection2(FFM, a, x, b, tau);
    }
    else {
        if (c - b > b - a) return GoldenSection2(FFM, a, b, x, tau);
        else return GoldenSection2(FFM, x, b, c, tau);
    }
    
}

double CGMethod::BinarySearch(ForceFieldManager& FFM, double a, double b )
{
    //double phi = 0.5 * (1 + sqrt(5)) ;
    double  eps = 0.001;
    
    while (fabs(b - a) > eps){
        
        double half_x1 = ((a + b)/2 - eps/4);
        double half_x2 = ((a + b)/2 + eps/4);
        if (FFM.ComputeEnergy(half_x1) <= FFM.ComputeEnergy(half_x2)) b = half_x2;
        else a = half_x1;
    }
     return (a + b) / 2;
}


