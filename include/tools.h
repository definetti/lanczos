#ifndef TOOLS_H
#define TOOLS_H

#include <stdlib.h>
#include <armadillo>
#include <time.h>


//Parameters of the system

using namespace arma;


    extern int n_sites;
    extern int n_electrons;
    extern int n_iter;
    extern double t;
    extern double U;

struct couple                       //a type consisting in a couple of variables, which will be useful
{
	int up;
	int down;
};

class Tools
{
    public:
    couple spins;
    ivec vector1; //service vector
    ivec vector2; //service vector
    imat H; //non diagonal hamiltonian
    vec prev; //old lanczos vector
    vec curr; //current lanczos vector
    vec next; //new lanczos vector
    mat lancmat; //lanczos matrix

    int binomial(int n, int k);
    int factorial (int n);

    ivec dec_to_bin(int a, ivec v);
    int filling (ivec v);
    int count_double(ivec v, ivec w);
    int are_equal(ivec v, ivec w);
    int are_neighbours(ivec v, ivec w);

    couple get_addresses(int i);
    int retrieve_addresses(couple v);

    void generate_H();
    void generate_next();
    vec build_lancmat();

};
#endif // TOOLS_H
