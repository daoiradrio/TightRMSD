#ifndef MOLECULE_HPP
#define MOLECULE_HPP

#include <atom.hpp>
#include <bond.hpp>
#include <utils.hpp>

#include <vector>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <Eigen/Dense>

using Atom_ptr = std::shared_ptr<Atom>;
using Bond_ptr = std::shared_ptr<Bond>;



class Molecule{

    public:
        Molecule();
        ~Molecule();

        int                     n_atoms;
        std::vector<Atom_ptr>   atoms;
        std::vector<Bond_ptr>   bonds;
        Eigen::MatrixX3d        coords;

        int                     compute_structure(std::string filepath);
        int                     read_xyz(std::string filepath);
    
    private:
        void                    set_connectivity();
        int                     set_bond_order(const Atom_ptr& atom1, const Atom_ptr& atom2);
};



#endif
