#ifndef ATOM_HPP
#define ATOM_HPP

#include <string>
#include <vector>

#define DIMENSION 3



struct Atom{
    std::string         element;
    int                 pse_num;
    int                 index;
    std::vector<double> coords;
    std::vector<int>    bond_partners;

    Atom(): coords(DIMENSION, 0.0) {}
    ~Atom(){}
};



#endif
