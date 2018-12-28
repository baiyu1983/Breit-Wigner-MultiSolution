#ifndef DRAW_H
#define DRAW_H
#include <vector>
#include "TComplex.h"
struct resonance{
    double mass;
    double width;
    double phase;
    double branchratio;
};

double Phi23(double sqrts,double M1,double M2,double M3,double M4,double M5);
double Phi23PiPiJPsi(double sqrts);

double amp2TF(double*, double *);
#endif
