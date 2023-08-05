#ifndef ATOM_HPP
#define ATOM_HPP

#include <string>
#include <vector>
#include <memory>

#define DIMENSION 3

struct Atom;
using Atom_ptr = std::shared_ptr<Atom>;



struct Atom{
    std::string             element;
    int                     pse_num;
    int                     index;
    std::vector<double>     coords;
    std::vector<Atom_ptr>   bond_partners;

    Atom(): coords(DIMENSION, 0.0) {}
    ~Atom(){}
};



#endif
