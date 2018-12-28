#ifndef DRAW_H
#define DRAW_H
#include <vector>
#include "TComplex.h"
/*
struct resonance{
    double mass;
    double width;
    double phase;
    double branchratio;
};
*/
//double Phi23(double sqrts,double M1,double M2,double M3,double M4,double M5);
//double Phi23PiPiJPsi(double sqrts);

double amp2TF_BW(double*, double *);
double amp2TF_BKG1(double*, double *);
double amp2TF_BKG2(double*, double *);
double amp2TF_Comb1(double *,double *);
double amp2TF_Comb2(double *,double *);
#endif
