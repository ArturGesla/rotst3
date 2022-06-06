#include <iostream>
#include "rs.hpp"
#include "rsStab.hpp"
#include "bodewadtWithPressure.hpp"
#include "bodewadtWithPressureStab.hpp"


int main(int argc, char *argv[])
{
    // rs1D::rs1D(51,-1);
    //  rs1DStab::rs1DStab();

    // if (argc == 2)
    // {
        
    //         std::cout << "void rs1DStab(double Re = 1.0, int n = 101, double Ro =1, double zmax = 30, double sigmaR=1, double sigmaI=0,double alphaR=0, double alphaI=0, double betaR=0, double betaI=0, int nev=10);"<<std::endl;
    //         std::cout << "./rotst 1.0 101 1 30 1.0 0.0 0.0 0.0 0.0 0.0 10 " << std::endl;
        
    // }
    // else 
    // {
        
    // std::cout << "Args: " << argc << std::endl;
    // for (int i = 1; i < argc; i++)
    // {
    //     std::cout << std::stof(argv[i]) << "\t";
    // }
    // std::cout << std::endl;
    // }
    
    // if (argc == 12)
    // {

    //     rs1DStab::rs1DStab(std::stof(argv[1]), std::stof(argv[2]), std::stof(argv[3]), std::stof(argv[4]),
    //                        std::stof(argv[5]), std::stof(argv[6]), std::stof(argv[7]), std::stof(argv[8]), std::stof(argv[9]),
    //                        std::stof(argv[10]), std::stof(argv[11]));
    // }
    // else if (argc !=2)
    // {
    //     rs1DStab::rs1DStab();
    // }

    // bwp::bwpBaseFlow();
    bwp::bwpStab();

    // rs1DStab::rs1DStab();

    return 0;
}