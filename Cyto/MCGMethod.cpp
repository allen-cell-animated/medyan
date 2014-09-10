//
//  MCGMethod.cpp
//  Cyto
//
//  Created by James Komianos on 9/10/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#include "MCGMethod.h"
using namespace std;


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
        
        (*it).coordinate[0] = (*it).coordinate[0] - d* (*it).force[0];
        (*it).coordinate[1] = (*it).coordinate[1] - d* (*it).force[1];
        (*it).coordinate[2] = (*it).coordinate[2] - d* (*it).force[2];
        
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
}



double CGMethod::GoldenSection(ForceFieldManager& FFM)
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
