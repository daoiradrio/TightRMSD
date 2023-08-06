#include <molecule.hpp>



Molecule::Molecule(){}



Molecule::~Molecule(){}



int Molecule::compute_structure(std::string filepath){
    int res;
    res = this->read_xyz(filepath);
    if (!res){
        return 0;
    }
    this->set_connectivity();
    return 1;
}



int Molecule::read_xyz(std::string filepath)
{
    std::shared_ptr<Atom> new_atom;
    std::string element;
    std::string new_label;
    std::string line;
    std::ifstream file;
    double energy;
    double xcoord, ycoord, zcoord;
    int atom_index;
    int line_index;

    this->n_atoms = 0;
    this->atoms.clear();
    this->coords.setZero();

    file.open(filepath);

    if (file.is_open()){
        line_index = 0;
        atom_index = 0;
        while (getline(file, line)){
            std::stringstream linestream(line);
            if (line_index == 0){
                linestream >> this->n_atoms;
                this->coords.resize(this->n_atoms, 3);
            }
            //else if (line_index == 1 && read_energy){
            //    linestream >> dummy1 >> dummy2 >> energy;
            //    this->energy = energy;
            //}
            else if (line_index >= 2){
                linestream >> element >> xcoord >> ycoord >> zcoord;  
                new_atom = std::make_shared<Atom>();
                new_atom->element = element;
                new_atom->pse_num = element_numbers[element];
                new_atom->index = atom_index;
                new_atom->coords = {xcoord, ycoord, zcoord};
                this->atoms.push_back(std::move(new_atom));

                this->coords(atom_index, 0) = xcoord;
                this->coords(atom_index, 1) = ycoord;
                this->coords(atom_index, 2) = zcoord;
                
                atom_index++;
            }
            line_index++;
        }
    }
    else{
        std::cout << "FAILED OPENING .xyz FILE!" << std::endl;
        return 0;
    }
    file.close();

    return 1;
}



void Molecule::set_connectivity()
{
    int                             bond_order;
    Atom_ptr                        atom1;
    Atom_ptr                        atom2;
    Bond_ptr                        new_bond;
    std::vector<Atom_ptr>::iterator i;
    std::vector<Atom_ptr>::iterator j;

    this->bonds.clear();

    for (i = this->atoms.begin(); i != this->atoms.end(); i++){
        for (j = std::next(i); j != this->atoms.end(); j++){
            atom1 = (*i);
            atom2 = (*j);
            bond_order = this->set_bond_order(atom1, atom2);    
            if (bond_order){
                atom1->bond_partners.push_back(atom2->index);
                atom2->bond_partners.push_back(atom1->index);
                new_bond = std::make_shared<Bond>();
                new_bond->atom1 = atom1;
                new_bond->atom2 = atom2;
                new_bond->bond_order = bond_order;
                this->bonds.push_back(new_bond);
            }
        }
    }

    /*
    for (Atom_ptr atom1: this->atoms){
        for (Atom_ptr atom2: this->atoms){
            bond_order = this->set_bond_order(atom1, atom2);    
            if (bond_order){
                atom1->bond_partners.push_back(atom2->index);
                atom2->bond_partners.push_back(atom1->index);
                new_bond = std::make_shared<Bond>();
                new_bond->atom1 = atom1;
                new_bond->atom2 = atom2;
                new_bond->bond_order = bond_order;
                this->bonds.push_back(new_bond);
            }
        }
    }
    */

    return;
}



int Molecule::set_bond_order(const Atom_ptr& atom1, const Atom_ptr& atom2)
{
    int bond_order = 0;
    double tolerance = 0.08;
    double single_bond = -1000.0;
    double double_bond = -1000.0;
    double triple_bond = -1000.0;

    double valence_radius_single1 = valence_radii_single[atom1->element];
    double valence_radius_single2 = valence_radii_single[atom2->element];

    double valence_radius_double1 = valence_radii_double[atom1->element];
    double valence_radius_double2 = valence_radii_double[atom2->element];

    double valence_radius_triple1 = valence_radii_triple[atom1->element];
    double valence_radius_triple2 = valence_radii_triple[atom2->element];

    double distance = sqrt(
        pow((atom1->coords[0] - atom2->coords[0]), 2) +
        pow((atom1->coords[1] - atom2->coords[1]), 2) +
        pow((atom1->coords[2] - atom2->coords[2]), 2)
    );

    if (valence_radius_single1 && valence_radius_single2){
        single_bond = valence_radius_single1 + valence_radius_single2 + tolerance;
    }
    if (valence_radius_double1 && valence_radius_double2){
        double_bond = valence_radius_double1 + valence_radius_double2 + tolerance;
    }
    if (valence_radius_triple1 && valence_radius_triple2){
        triple_bond = valence_radius_triple1 + valence_radius_triple2 + tolerance;
    }

    if (distance <= triple_bond){
        bond_order = 3;
    }
    else if (distance <= double_bond){
        bond_order = 2;
    }
    else if (distance <= single_bond){
        bond_order = 1;
    }

    return bond_order;
}
