#include <iostream>
#include <Eigen/Dense>
#include <utility.hpp>



int main(int argc, char* argv[])
{
    if (argc != 3){
        std::cout << "\nINVALID INPUT!\n" << std::endl;
        return 1;
    }
    std::string file1 = argv[1];
    std::string file2 = argv[2];

    Eigen::MatrixX3d matched_coords1;
    Eigen::MatrixX3d matched_coords2;

    std::vector<Atom_ptr> atoms1 = compute_structure(file1);
    std::vector<Atom_ptr> atoms2 = compute_structure(file2);

    match_atoms(atoms1, atoms2, matched_coords1, matched_coords2);

    kabsch(matched_coords1, matched_coords2);

    std::cout << rmsd(matched_coords1, matched_coords2) << std::endl;

    //std::cout << matched_coords1.rows() << " " << matched_coords1.cols() << std::endl;
    //std::cout << matched_coords2.rows() << " " << matched_coords2.cols() << std::endl;
    
    return 0;
}
