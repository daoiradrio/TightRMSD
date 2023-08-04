#include <molecule.hpp>



Molecule::Molecule(){}


Molecule::~Molecule(){}


void Molecule::compute_structure(std::string filepath){
    this->read_xyz(filepath);
    this->set_connectivity();
    return;
}


int Molecule::read_xyz(std::string filepath)
{
    return 0;
}


void Molecule::set_connectivity()
{
    return;
}


int Molecule::set_bond_order(int i, int j)
{
    return 0;
}
