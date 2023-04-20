#include "build.h"
#include "query.h"
#include <iostream>
#include <fstream>
// #include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Dense>
using namespace std;

int main()
{
    Eigen::MatrixXf m = Eigen::MatrixXf::Random(3,2);
     std::cout << "Here is the matrix m:" << std::endl << m << std::endl;
    Eigen::JacobiSVD<Eigen::MatrixXf, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(m);
     cout << "Its singular values are:" << endl << svd.singularValues() << endl;
     cout << "Its left singular vectors are the columns of the thin U matrix:" << endl << svd.matrixU() << endl;
     cout << "Its right singular vectors are the columns of the thin V matrix:" << endl << svd.matrixV() << endl;
    Eigen::Vector3f rhs(1, 0, 0);
     cout << "Now consider this rhs vector:" << endl << rhs << endl;
     cout << "A least-squares solution of m*x = rhs is:" << endl << svd.solve(rhs) << endl;
    Solid *s = torus(1, 100, 10, 100, 100);
    // Solid *s = circle(1, 0.0, 100.0, 10, 0.0, 10);
    // Solid *s = block(1, 10, 10, 10);
    triangel t = delaunay(s);

    std::cout << "-----------" << std::endl;

    //write doc
    std::string dir("view.obj");
    // std::ofstream outfile(dir);
    std::ofstream outfile(dir);
    

    //wirte vetex
    auto it = t.vtxarry.begin();
    auto end = t.vtxarry.end();
    for (; it != end; ++it)
    {
        outfile << "v " << (*it)->vcoord[0] << " " <<
            (*it)->vcoord[1] << " " << (*it)->vcoord[2] << std::endl;
    }

    // write face
    auto itf = t.face.begin();
    auto endf = t.face.end();
    for (; itf != endf; ++itf)
    {
        outfile << "f";
        std::vector<int> facevecter = *itf;
        for (auto item : facevecter)
        {
            outfile << " " << item;
        }
        outfile << std::endl;
    }

    outfile.close();
    
    // listsolid(s);
    // std::cout << 1 <<std::endl;
    return 0;
}