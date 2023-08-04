#ifndef MOLECULE_HPP
#define MOLECULE_HPP

#include <atom.hpp>
#include <bond.hpp>

#include <vector>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <Eigen/Dense>



class Molecule{

    public:
        Molecule();
        ~Molecule();

        int                                 n_atoms;
        double                              energy = 0;
        std::vector<std::shared_ptr<Atom>>  atoms;
        std::vector<std::shared_ptr<Bond>>  bonds;
        Eigen::MatrixX3d                    coords;

        void    compute_structure(std::string filepath);
        int     read_xyz(std::string filepath);
    
    private:
        void    set_connectivity();
        int     set_bond_order(int i, int j);
};


#endif
