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
const double r = 0.61803399;
const double c = 1 - r;

inline void CGMethod::swap(double &a, double &b) {
    double tmp = a;
    a = b;
    b = tmp;
}

inline void CGMethod::shift2(double &a, double &b, double c) {
    a = b;
    b = c;
}

inline void CGMethod::shift3(double &a, double &b, double &c, double d){
    a = b;
    b = c;
    c = d;
}

inline double CGMethod::sign(double a, double b) {
    
    return b >= 0 ? fabs(a) : -fabs(a);
    
}

void CGMethod::makeBracket(ForceFieldManager &FFM, double &ax, double &bx, double &cx, double &fa, double &fb, double &fc) {
    
    const double GLIMIT = 100.0, TINY = 1.0e-20;
    double ulim, u, r, q, fu;
    
    fa = FFM.ComputeEnergy(ax);
    fb = FFM.ComputeEnergy(bx);
    
    ///go in correct direction
    if (fb > fa){
        swap(ax, bx);
        swap(fb, fa);
    }

    //first guess
    cx = bx + phi * (bx - ax);
    fc = FFM.ComputeEnergy(cx);
    
    while(fb > fc) {
        
        r = (bx - ax) * (fb - fc);
        q = (bx - cx) * (fb - fa);
        u = bx - ((bx - cx) * q - (bx - ax) * r) / (2.0 * sign(max(fabs(q - r), TINY), q - r));
        
        ulim = bx + GLIMIT * (cx - bx);
        
        if ((bx - u) * (u - cx) > 0.0) {
            fu = FFM.ComputeEnergy(u);
            if (fu < fc) {
                ax = bx;
                bx = u;
                fa = fb;
                fb = fu;
                return;
            }
            else if (fu > fb) {
                cx = u;
                fc = fu;
                return;
            }
            u = cx + phi * (cx - bx);
            fu = FFM.ComputeEnergy(u);
        }
        else if ((cx - u) * (u - ulim) > 0.0) {
            
            fu = FFM.ComputeEnergy(u);
            if (fu < fc) {
                shift3(bx, cx, u, cx + phi * (cx - bx));
                shift3(fb, fc, fu, FFM.ComputeEnergy(u));
            }
        }
        else if((u - ulim) * (ulim - cx) >= 0.0) {
            
            u = ulim;
            fu = FFM.ComputeEnergy(u);
        }
        else {
            u = cx + phi * (cx - bx);
            fu = FFM.ComputeEnergy(u);
        }
        shift3(ax, bx, cx, u);
        shift3(fa, fb, fc, fu);
    }
}


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
	double b = 15;
    double phi = 0.5 * (1 + sqrt(5) );
    double inv_phi = 1/phi;
	double x1 = b - inv_phi * (b - a);
	double x2 = a + inv_phi * (b - a);
    
	while (fabs(b - a) > EPS)
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

double CGMethod::GoldenSection2(ForceFieldManager& FFM, double ax, double bx, double cx, double tau)
{
    double x;
    if (cx - bx > bx - ax)
        x = bx + resphi * (cx - bx);
    else
        x = bx - resphi * (bx - ax);
    if (fabs(cx - ax) < tau * (fabs(bx) + fabs(x)))
        return (cx + ax) / 2;
    
    double fx = FFM.ComputeEnergy(x);
    double fb = FFM.ComputeEnergy(bx);
    
    if (fx < fb) {
        if (cx - bx > bx - ax) return GoldenSection2(FFM, bx, x, cx, tau);
        else return GoldenSection2(FFM, ax, x, bx, tau);
    }
    else {
        if (cx - bx > bx - ax) return GoldenSection2(FFM, ax, bx, x, tau);
        else return GoldenSection2(FFM, x, bx, cx, tau);
    }
    
}

double CGMethod::GoldenSection3(ForceFieldManager& FFM, double ax, double bx, double cx, double tol) {
    
    double f1, f2, x0, x1, x2, x3;
    
    x0 = ax;
    x3 = cx;
    if (fabs(cx - bx) > fabs(bx - ax)) {
        x1 = bx;
        x2 = bx + c * (cx - bx);
    }
    else {
        x2 = bx;
        x1 = bx - c * (bx - ax);
    }
    f1 = FFM.ComputeEnergy(x1);
    f2 = FFM.ComputeEnergy(x2);
    
    while (fabs(x3 - x0) > tol * (fabs(x1) + fabs(x2))) {
        if (f2 < f1) {
            shift3(x0, x1, x2, r * x2 + c * x3);
            shift2(f1, f2, FFM.ComputeEnergy(x2));
        }
        else {
            shift3(x3, x2, x1, r * x1 + c * x0);
            shift2(f2, f1, FFM.ComputeEnergy(x1));
        }
    }
    
    if(f1 < f2) return x1;
    else return x2;
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


