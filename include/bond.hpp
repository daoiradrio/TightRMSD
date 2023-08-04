#ifndef BOND_HPP
#define BOND_HPP

#include <atom.hpp>

#include <memory>



struct Bond{
    std::shared_ptr<Atom>   atom1;
    std::shared_ptr<Atom>   atom2;
    int                     bond_order = 0;
};



#endif
