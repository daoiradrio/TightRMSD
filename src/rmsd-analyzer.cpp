#include <rmsd-analyzer.hpp>



RMSDAnalyzer::RMSDAnalyzer(){};



RMSDAnalyzer::~RMSDAnalyzer(){};



void RMSDAnalyzer::compute_tight_rmsd(std::string file1, std::string file2)
{
    return;
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
