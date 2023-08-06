#include <iostream>
#include <molecule.hpp>
#include <rmsd-analyzer.hpp>



int main(int argc, char* argv[]){
    //Molecule mol1;
    //Molecule mol2;
    RMSDAnalyzer analyzer;
    std::string file1;
    std::string file2;

    if (argc != 3){
        std::cout << "\nINVALID INPUT!\n" << std::endl;
        return 1;
    }
    file1 = argv[1];
    file2 = argv[2];

    analyzer.compute_tight_rmsd(file1, file2);
    
    return 0;
}
