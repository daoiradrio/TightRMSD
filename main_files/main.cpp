#include <iostream>
#include <molecule.hpp>
#include <rmsd-analyzer.hpp>



int main(){
    Molecule mol;
    RMSDAnalyzer analyzer;

    //mol.read_xyz("/home/dario/TightRMSD/inputfiles/Tyrosin.xyz");
    mol.compute_structure("/home/dario/TightRMSD/inputfiles/Tyrosin.xyz");
    
    return 0;
}
