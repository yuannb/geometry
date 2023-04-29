#include "build.h"
#include "query.h"
#include <iostream>
#include <fstream>
#include "planeequ.h"
// #include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "split.h"
using namespace std;

int main()
{
    std::shared_ptr<Solid> s = torus(1, 100, 10, 100, 100);
    std::shared_ptr<Solid> above = nullptr;
    std::shared_ptr<Solid> below = nullptr;
    Eigen::Vector4d sp{ 1, 0, 0, 0 };
    // std::shared_ptr<Solid> s = block(1, 10, 10, 10);
    split(s, sp, above, below);
    // // Solid *s = circle(1, 0.0, 100.0, 10, 0.0, 10);
    triangel t = delaunay(above);

    std::cout << "-----------" << std::endl;

    // write doc
    std::string dir("view.obj");
    // std::ofstream outfile(dir);
    std::ofstream outfile(dir);


    // wirte vetex
    auto it = t.vtxarry.begin();
    auto end = t.vtxarry.end();
    for (; it != end; ++it)
    {
       outfile << "v " << (*it)->vcoord[0] << " " <<
           (*it)->vcoord[1] << " " << (*it)->vcoord[2] << std::endl;
    }

    //write face
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

    // write doc
    dir = std::string("view1.obj");
    // std::ofstream outfile(dir);
    std::ofstream outfile2(dir);


    // wirte vetex
    triangel t2 = delaunay(below);
    it = t2.vtxarry.begin();
    end = t2.vtxarry.end();
    for (; it != end; ++it)
    {
       outfile2 << "v " << (*it)->vcoord[0] << " " <<
           (*it)->vcoord[1] << " " << (*it)->vcoord[2] << std::endl;
    }

    // write face
    itf = t2.face.begin();
    endf = t2.face.end();
    for (; itf != endf; ++itf)
    {
       outfile2 << "f";
       std::vector<int> facevecter = *itf;
       for (auto item : facevecter)
       {
           outfile2 << " " << item;
       }
       outfile2 << std::endl;
    }

    outfile2.close();

     listsolid(s);
     std::cout << 1 <<std::endl;
    return 0;
}