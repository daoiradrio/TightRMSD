#ifndef RMSDANALYZER_HPP
#define RMSDANALYZER_HPP

#include <molecule.hpp>

#include <string>
#include <vector>
#include <queue>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/SVD>

#include <iostream>



class RMSDAnalyzer{

    public:

        RMSDAnalyzer();
        ~RMSDAnalyzer();

        void    compute_tight_rmsd(std::string file1, std::string file2);

        std::vector<std::vector<int>>   spheres(Molecule& molecule, Atom_ptr start_atom);

    private:

        double                          rmsd(Eigen::MatrixX3d coords1, Eigen::MatrixX3d coords2);
        void                            kabsch(Eigen::MatrixX3d& coords1, Eigen::MatrixX3d& coords2);
        void                            match_atoms(std::vector<Atom_ptr> atoms1, std::vector<Atom_ptr> atoms2);
        //std::vector<std::vector<int>>   spheres(Molecule& molecule, Atom_ptr start_atom);
};



#endif
