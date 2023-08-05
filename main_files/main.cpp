#include <iostream>
#include <molecule.hpp>
#include <rmsd-analyzer.hpp>



int main(){
    Molecule mol;

    mol.compute_structure("../inputfiles/Tyrosin.xyz");
    
    return 0;
}
