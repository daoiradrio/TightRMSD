#ifndef BOND_HPP
#define BOND_HPP

#include <atom.hpp>

#include <memory>

using Atom_ptr = std::shared_ptr<Atom>;



struct Bond{
    Atom_ptr    atom1;
    Atom_ptr    atom2;
    int         bond_order = 0;
};



#endif
