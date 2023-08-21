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

    std::cout << std::endl;
    std::cout << "File 1: " << file1 << std::endl;
    std::cout << "File 2: " << file2 << std::endl;
    std::cout << std::endl;
    std::cout << "RMSD: " << tight_rmsd(file1, file2) << std::endl;
    std::cout << std::endl;
    
    return 0;
}
