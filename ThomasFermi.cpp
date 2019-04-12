/* C++ Program by Josh Frechem
   Compile with g++ with OpenMP support or
   On MacOS:  clang++ -Xpreprocessor -fopenmp -lomp -o Scattering ThomasFermi.cpp    */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <complex>
#include <string>
#include <cstring>
#include <string.h>
#include <cstdlib>
#include <omp.h>
#include <chrono>
#include <ctime>

#include "Headers/Globals.h"
#include "Headers/Besselj.h"
#include "Headers/Numerov.h"
#include "Headers/Potentials.h"


//Units are eV,m,s,
//---------------------------------------------------------------------------------
//Physical Constants---------------------------------------------------------------  
const double PI = 3.141592653589793238462643383279502884;
const double h2m = 7.6359;                  //eV-A^2
const double m=0.5109989461*pow(10,6);      //electron mass eV/c^2
const double ee=14.409;                     //electron charge squared in eV-A
const double aBohr=h2m/ee;                  //Bohr radius in A
const double hbar=6.582119514*pow(10,-16);  //hbar eV-s
const double ch=1.6021766208*pow(10,-19);   //electron charge in C
const double SOL=299792458;                 //Speed of light in m/s
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------


// Begin Simulation main function
int main()
{
    int MeshSize=2000;  //mesh size
    int nMax=50;        //How far to calculate Bessels out to
    int Z=10;           //atomic number
    int Energy=5000;    //energy level in eV to calculate (number of nodes for numerov)
    double Rmax=2.0*pow(10,-10); //Mesh radius limit in angstroms
    double h=Rmax/MeshSize;      //step size in angstroms
    double k=sqrt(2*m*Energy/hbar/hbar/SOL/SOL);                 //wavenumber 1/m
    double cval=-2*m*(2*Rmax)*(2*Rmax)*Energy/hbar/hbar/SOL/SOL; //constant for quadratic solving
    int lmax=(-1+sqrt(1-4*cval))/(2);                            //quadratically solving for lmax using E=l(l+1)*hbar^2/2mr^2


    //Timekeeping------------------------------------------------------------------------------------------------------------------------------------------------------
    std::chrono::time_point<std::chrono::system_clock> startTTM, endTTM, startDIFF,endDIFF, startTOT, endTOT;
    startTOT=std::chrono::system_clock::now();
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------


    //Setting up parallel processing-----------------------------------------------------------------------------------------------------------------------------------
        int numcores= omp_get_num_procs(); //set number of parallel workers by finding number of CPU workers. id1 and id2 are optimized for the amount of work per task
        int threads=1.0*(numcores);
        std::cout<<"\n           Thomas-Fermi Atomic Potential Scattering\n Shared Memory Parallelism Enabled:\t"<<numcores<<" processing cores used"<<std::endl;
        std::cout<<"_________________________________________________________________________________________________________"<<std::endl;
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------

    std::vector<std::vector<double > > BesseljSol (MeshSize+1,std::vector<double>(nMax+1));
    std::vector<std::vector<double > > BesselySol (MeshSize+1,std::vector<double>(nMax+1));
    std::vector<std::vector<double > > Potential (MeshSize+1,std::vector<double>(2));


    //=================================================================================================================================================================
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------
    //   This is the end of simulation setup and beginning of driving function
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------
    //=================================================================================================================================================================


            for (int i=0;i<=MeshSize;++i)
            {
                BesseljSol[i][0]=Rmax/MeshSize*i;  //First column is argument array
                BesselySol[i][0]=Rmax/MeshSize*i;  //First column is argument array
                Potential[i][0]=Rmax/MeshSize*i;  //First column is argument array
                Potential[i][1]=ThomasFermiPotential(Z,Potential[i][0]); //Second column is V array
                //Potential[i][1]=HarmonicPotential(Potential[i][0],1.0,1.0);
                //Potential[i][1]=ZeroPotential(Z,Potential[i][0]);
            }

            omp_set_num_threads(threads);
            #pragma omp parallel for shared(nMax, BesseljSol) //GenerateBessel Functions Recursively in Parallel
            for(int i=0;i<=MeshSize;i++)
            {
                for (int n=1; n<nMax+1; ++n)
                {
                    BesseljSol[i][n]=BesselJdown(n-1,k*BesseljSol[i][0]);
                    BesselySol[i][n]=BesselJdown(n-1,k*BesselySol[i][0]);
                }
            }

            //double pp=numerov(Energy,MeshSize,h,Potential);
            //std::cout<<pp<<std::endl;


    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------
    //    This is the end of the driving function and beginning of writing out data to file and/or console
    //------------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------

    std::ofstream Solutions_out ("BesseljSol.txt");
      if (Solutions_out.is_open())
      {
            for (int b=0; b<=MeshSize; ++b)
            {
                for (int c=0; c<nMax+1; ++c)
                {
                    Solutions_out<<BesseljSol[b][c]<<"\t";
                }
                Solutions_out<<Potential[b][1]<<"\n";
            }
            Solutions_out.close();
      }





        endTOT=std::chrono::system_clock::now();
        std::chrono::duration<double> TOTtime=endTOT-startTOT;
        std::cout<<"\nTotal execution time: "<<TOTtime.count()<<std::endl<<std::endl<<"End of line, man."<<std::endl<<std::endl<<std::endl;

    std::vector<std::vector<double > >().swap(BesseljSol);
    std::vector<std::vector<double > >().swap(Potential);

    return 0;
}



//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//    End of line, man
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
