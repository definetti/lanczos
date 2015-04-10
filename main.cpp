#include <iostream>
#include <iomanip>      // std::setprecision
#include <fstream>
#include <stdlib.h>
#include <armadillo>
#include <time.h>
#include "tools.h"

using namespace std;
using namespace arma;

int n_sites=10;
int n_iter=20;
double t=1.0;
double U=2.0;
int n_electrons;


int main(int argc, char *argv[])
{
    vec energies=zeros<vec>(n_iter);
    Tools stuff;

    n_electrons = n_sites/2;
    stuff.generate_H();
    energies=stuff.build_lancmat();

    cout.precision(10);
    cout<<"GS energy = "<<energies(0)<<endl;
    cout<<"Per site = "<<energies(0)/n_sites<<endl;

    return 0;
}
