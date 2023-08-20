#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <string>
#include <memory>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <queue>
#include <algorithm>
#include <utility>
#include <Eigen/Dense>
#include <Eigen/SVD>

#define DIMENSION 3

struct Atom;
using Atom_ptr = std::shared_ptr<Atom>;

std::vector<Atom_ptr>           compute_structure(std::string filepath);
bool                            connected(const Atom_ptr& atom1, const Atom_ptr& atom2);
std::vector<std::vector<int>>   spheres(const std::vector<Atom_ptr>& atoms, const Atom_ptr& start_atom);
void                            match_atoms(std::vector<Atom_ptr> atoms1, std::vector<Atom_ptr> atoms2);
void                            kabsch(Eigen::MatrixX3d& coords1, Eigen::MatrixX3d& coords2);



std::unordered_map<std::string, int> element_numbers(
    {
        {"H", 1},
        {"B", 5},
        {"C", 6},
        {"N", 7},
        {"O", 8},
        {"F", 9},
        {"Cl", 17},
        {"Br", 35},
        {"I", 53}
    }
);



std::unordered_map<std::string, double> valence_radii_single(
    {
        {"C", 0.75},
        {"N", 0.71},
        {"O", 0.66},
        {"H", 0.32},
        {"B", 0.85},
        {"F", 0.64},
        {"Cl", 0.99},
        {"Br", 1.14},
        {"I", 1.33}
    }
);



struct Atom
{
    std::string             element;
    int                     pse_num;
    int                     index;
    Eigen::RowVector3d      coords;
    std::vector<int>        bond_partners;
    std::vector<Atom_ptr>   eq_atoms;
};



std::vector<Atom_ptr> compute_structure(std::string filepath)
{
    std::vector<Atom_ptr>   atoms;
    std::shared_ptr<Atom>   new_atom;
    std::string             line;
    std::ifstream           file;
    std::string             element;
    double                  xcoord, ycoord, zcoord;
    int                     atom_index;
    int                     line_index;

    file.open(filepath);

    if (file.is_open()){
        line_index = 0;
        atom_index = 0;
        while (getline(file, line)){
            if (line_index >= 2){
                std::stringstream linestream(line);
                linestream >> element >> xcoord >> ycoord >> zcoord;  
                new_atom = std::make_shared<Atom>();
                new_atom->element = element;
                new_atom->pse_num = element_numbers[element];
                new_atom->index = atom_index;
                new_atom->coords << xcoord, ycoord, zcoord;
                atoms.push_back(std::move(new_atom));
                atom_index++;
            }
            line_index++;
        }
    }
    else{
        std::cout << "\nFAILED OPENING .xyz FILE!\n" << std::endl;
        return {};
    }
    file.close();

    Atom_ptr                        atom1;
    Atom_ptr                        atom2;
    std::vector<Atom_ptr>::iterator i;
    std::vector<Atom_ptr>::iterator j;

    for (i = atoms.begin(); i != atoms.end(); i++){
        for (j = std::next(i); j != atoms.end(); j++){
            atom1 = (*i);
            atom2 = (*j);   
            if (connected(atom1, atom2)){
                atom1->bond_partners.push_back(atom2->index);
                atom2->bond_partners.push_back(atom1->index);
            }
        }
    }

    return atoms;
}



bool connected(const Atom_ptr& atom1, const Atom_ptr& atom2)
{
    double tolerance = 0.08;
    double single_bond = -1000.0;

    double valence_radius_single1 = valence_radii_single[atom1->element];
    double valence_radius_single2 = valence_radii_single[atom2->element];

    double distance = sqrt(
        pow((atom1->coords(0) - atom2->coords(0)), 2) +
        pow((atom1->coords(1) - atom2->coords(1)), 2) +
        pow((atom1->coords(2) - atom2->coords(2)), 2)
    );

    if (valence_radius_single1 && valence_radius_single2){
        single_bond = valence_radius_single1 + valence_radius_single2 + tolerance;
    }

    if (distance <= single_bond){
        return true;
    }

    return false;
}



// WIRD WAHRSCHEINLICH RETURN VALUE BRAUCHEN
void match_atoms(
    const std::vector<Atom_ptr>&    atoms1,
    const std::vector<Atom_ptr>&    atoms2,
    Eigen::MatrixX3d&               coords1,
    Eigen::MatrixX3d&               coords2
){   
    int                                         i, j;
    double                                      deviation;
    double                                      min;
    Atom_ptr                                    min_atom;
    Eigen::VectorXd                             d_ref;
    Eigen::VectorXd                             d;

    i = 0;
    coords1.resize(atoms1.size(), 3);
    coords2.resize(atoms2.size(), 3);

    for (Atom_ptr atom1: atoms1){
        for (Atom_ptr atom2: atoms2){
            if (spheres(atoms1, atom1) == spheres(atoms2, atom2)){
                atom1->eq_atoms.push_back(atom2);
                atom2->eq_atoms.push_back(atom1);
            }
        }
        if (atom1->eq_atoms.size() == 1){
            coords1.block<1,3>(i,0) = atom1->coords;
            coords2.block<1,3>(i,0) = atom1->eq_atoms[0]->coords;
            i++;
        }
        else if (atom1->eq_atoms.empty()){
            // HIER SECURITY OPTION (REIHENFOLGE UNVERÄNDERT ZURÜCKGEBEN)
            continue;
        }
    }

    d_ref.resize(coords1.rows());
    d.resize(coords2.rows());
    
    for (Atom_ptr atom1: atoms1){
        if (atom1->eq_atoms.size() == 1){
            continue;
        }
        for (j = 0; j < coords1.rows(); j++){
            d_ref(j) = (atom1->coords - coords1.row(j)).norm();
        }
        min = 1000000.0;
        for (Atom_ptr atom2: atom1->eq_atoms){
            for (j = 0; j < coords2.rows(); j++){
                d(j) = (atom2->coords - coords2.row(j)).norm();
            }
            deviation = (d - d_ref).norm();
            if (deviation < min){
                min = deviation;
                min_atom = atom2;
            }
        }
        coords1.block<1,3>(i,0) = atom1->coords;
        coords2.block<1,3>(i,0) = min_atom->coords;
        i++;
    }

    return;
}



std::vector<std::vector<int>> spheres(const std::vector<Atom_ptr>& atoms, const Atom_ptr& start_atom)
{
    std::vector<std::vector<int>>   spheres;
    std::vector<int>                new_sphere;
    std::queue<Atom_ptr>            next_atoms;
    Atom_ptr                        next_atom;
    int                             current_sphere_iter;
    int                             next_sphere_iter;
    std::vector<int>                memory;

    current_sphere_iter = start_atom->bond_partners.size();
    next_sphere_iter = 0;
    memory = {start_atom->index};

    for (int neighbor_index: start_atom->bond_partners){
        next_atoms.push(atoms[neighbor_index]);
    }

    while (!next_atoms.empty()){
        new_sphere = {};
        while (current_sphere_iter){
            next_atom = next_atoms.front();
            next_atoms.pop();
            new_sphere.push_back(next_atom->pse_num);
            current_sphere_iter--;
            for (int neighbor_index: next_atom->bond_partners){
                if (std::find(memory.begin(), memory.end(), neighbor_index) != memory.end()){
                    continue;
                }
                memory.push_back(neighbor_index);
                next_atoms.push(atoms[neighbor_index]);
                next_sphere_iter++;
            }
        }
        std::sort(new_sphere.begin(), new_sphere.end());
        spheres.push_back(new_sphere);
        current_sphere_iter = next_sphere_iter;
        next_sphere_iter = 0;
    }

    return spheres;
}



void kabsch(Eigen::MatrixX3d& coords1, Eigen::MatrixX3d& coords2)
{
    int                 i;
    double              det;
    Eigen::Matrix3Xd    coords1_T;
    Eigen::Matrix3Xd    coords2_T;
    Eigen::MatrixX3d    H;
    Eigen::Matrix3d     helper_mat;
    Eigen::Matrix3d     R;
    Eigen::Vector3d     center1;
    Eigen::Vector3d     center2;
    
    center1 = {0.0, 0.0, 0.0};
    center2 = {0.0, 0.0, 0.0};
    for (i = 0; i < coords1.rows(); i++){
        center1 = center1 + coords1.row(i).transpose();
        center2 = center2 + coords2.row(i).transpose();
    }
    
    center1 = (1.0/(double)coords1.rows()) * center1;
    center2 = (1.0/(double)coords2.rows()) * center2;

    for (i = 0; i < coords1.rows(); i++){
        coords1.row(i) = coords1.row(i) - center1.transpose();
        coords2.row(i) = coords2.row(i) - center2.transpose();
    }
    
    coords1_T = coords1.transpose();

    H = coords1_T * coords2;

    //Eigen::JacobiSVD<Eigen::MatrixX3d> svd(H, Eigen::ComputeFullV | Eigen::ComputeFullU);
    Eigen::BDCSVD<Eigen::MatrixX3d> svd(H, Eigen::ComputeFullV | Eigen::ComputeFullU);

    det = (svd.matrixV() * svd.matrixU().transpose()).determinant();
    
    if (det >= 0.0){
        det = 1.0;
    }
    else{
        det = -1.0;
    }

    helper_mat << 1.0, 0.0, 0.0,
         0.0, 1.0, 0.0,
         0.0, 0.0, det;

    R = (svd.matrixV() * helper_mat) * svd.matrixU().transpose();

    for (i = 0; i < coords2.rows(); i++){
        coords2.row(i) = coords2.row(i) * R;
    }

    return;
}



#endif
