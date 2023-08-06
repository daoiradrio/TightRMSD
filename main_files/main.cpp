#include <iostream>
#include <utility.hpp>



int main(int argc, char* argv[])
{
    if (argc != 3){
        std::cout << "\nINVALID INPUT!\n" << std::endl;
        return 1;
    }
    std::string file1 = argv[1];
    std::string file2 = argv[2];

    std::vector<Atom_ptr> atoms1 = compute_structure(file1);
    std::vector<Atom_ptr> atoms2 = compute_structure(file2);
    
    return 0;
}
