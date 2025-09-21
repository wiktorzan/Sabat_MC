#ifndef ENERGY_SMEARER_H
#define ENERGY_SMEARER_H

class EnergySmearer {
    private:
        double a,b,c;
    public:
        EnergySmearer(double a_param, double b_param, double c_param);
        double SmearEnergy(double energy);
};

#endif