#include <rmsd-analyzer.hpp>



RMSDAnalyzer::RMSDAnalyzer(){};



RMSDAnalyzer::~RMSDAnalyzer(){};



void RMSDAnalyzer::compute_tight_rmsd(std::string file1, std::string file2)
{
    Molecule    mol1;
    Molecule    mol2;

    mol1.compute_structure(file1);
    mol2.compute_structure(file2);

    return;
}



void RMSDAnalyzer::match_atoms(std::vector<Atom_ptr> atoms1, std::vector<Atom_ptr> atoms2)
{
    std::vector<Atom_ptr>   eq_atoms;

    for (Atom_ptr atom1: atoms1){
        for (Atom_ptr atom2: atoms2){
            if (this->spheres(atoms1, atom1) == this->spheres(atoms2, atom2)){
                std::cout << "hier" << std::endl;
            }
        }
    }

    return;
}



std::vector<std::vector<int>> RMSDAnalyzer::spheres(const std::vector<Atom_ptr>& atoms, const Atom_ptr& start_atom)
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



double RMSDAnalyzer::rmsd(Eigen::MatrixX3d coords1, Eigen::MatrixX3d coords2)
{
    int             i;
    int             n = coords1.rows();
    double          rmsd = 0.0;
    Eigen::Vector3d diff_vec;

    for (i = 0; i < n; i++){
        diff_vec = coords1.row(i) - coords2.row(i);
        rmsd    += diff_vec.dot(diff_vec);
    }
    rmsd = (1.0/(double)n) * rmsd;
    rmsd = sqrt(rmsd);

    return rmsd;
}



void RMSDAnalyzer::kabsch(Eigen::MatrixX3d& coords1, Eigen::MatrixX3d& coords2)
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

    Eigen::JacobiSVD<Eigen::MatrixX3d> svd(H, Eigen::ComputeFullV | Eigen::ComputeFullU);
    //Eigen::BDCSVD<Eigen::MatrixX3d> svd(H, Eigen::ComputeFullV | Eigen::ComputeFullU);

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
