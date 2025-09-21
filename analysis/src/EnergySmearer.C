#include "EnergySmearer.h"
#include "TRandom3.h"

EnergySmearer::EnergySmearer(double a_param, double b_param, double c_param) 
    : a(a_param), b(b_param), c(c_param) {}

double EnergySmearer::SmearEnergy(double energy) {
    return gRandom->Gaus(0,1)*(a + b*sqrt(energy + c*pow(energy, 2)))/(2.35482004503);
}