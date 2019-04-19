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
#include "Headers/Bessely.h"
#include "Headers/LegendreP.h"
#include "Headers/Numerov.h"
#include "Headers/Potentials.h"
#include "Headers/Derivatives.h"


//Units are eV,m,s,
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

//Setting up parallel processing-----------------------------------------------------------------------------------------------------------------------------------
    int numcores= omp_get_num_procs(); //set number of parallel workers by finding number of CPU workers. id1 is to be optimized for the amount of work per task
    double id1=1.0;
    const int threads=id1*(numcores);
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------


// Begin Simulation main function
int main()
{
    int MeshSize=2000;  //mesh size
    int Z=18;           //atomic number
    int Energy=20;    //energy level in eV to calculate (number of nodes for numerov)
    double Rmax=2.0*pow(10,-10); //Mesh radius limit in angstroms
    double scale=5.0;            //Scale for how far out past Rmax to go
    double h=5*Rmax/MeshSize;    //step size in angstroms
    double k=sqrt(2*m*Energy/hbar/hbar/SOL/SOL);                 //wavenumber 1/m
    double cval=-2*m*(2*Rmax)*(2*Rmax)*Energy/hbar/hbar/SOL/SOL; //constant for quadratic solving
    int lmax=(-1+sqrt(1-4*cval))/(2);


    //Timekeeping------------------------------------------------------------------------------------------------------------------------------------------------------
        std::chrono::time_point<std::chrono::system_clock> startTTM, endTTM, startDIFF,endDIFF, startTOT, endTOT;
        startTOT=std::chrono::system_clock::now();
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------



        std::cout<<"\n           Lenz-Jensen Atomic Potential Scattering\n Shared Memory Parallelism Enabled:\t"<<numcores<<" processing cores used"<<std::endl;
        std::cout<<"_________________________________________________________________________________________________________"<<std::endl;
        std::cout<<"\nNumber of partial waves to use: "<<lmax<<std::endl<<std::endl;

    //Memory Allocation------------------------------------------------------------------------------------------------------------------------------------------------
        std::vector<std::vector<double > > BesseljSol (MeshSize+1,std::vector<double>(lmax+1,0));
        std::vector<std::vector<double > > BesselySol (MeshSize+1,std::vector<double>(lmax+1,0));
        std::vector<std::vector<double > > LegendreSol (181,std::vector<double>(lmax+1,0));
        std::vector<std::vector<double > > Potential (MeshSize+1,std::vector<double>(2,0));
        std::vector<double > Rarray (MeshSize,0.0);
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------


    //=================================================================================================================================================================
    //   This is the end of simulation setup and beginning of driving function
    //=================================================================================================================================================================


            for (int i=0;i<=MeshSize;++i)
            {
                Rarray[i]=scale*Rmax/MeshSize*i;  //First column is argument array
                Potential[i][0]=scale*Rmax/MeshSize*i;   //First column is argument array
                Potential[i][1]=ThomasFermiPotential(Z,Potential[i][0]); //Second column is V array
                //Potential[i][1]=HarmonicPotential(Potential[i][0],1.0,1.0);
                //Potential[i][1]=ZeroPotential(Z,Potential[i][0]);
            }

            omp_set_num_threads(threads);
            #pragma omp parallel for shared(MeshSize, lmax, BesseljSol, BesselySol) //Generate Bessel Functions Recursively in Parallel
            for(int i=0;i<=MeshSize;i++)
            {
                for (int n=0; n<=lmax; ++n)
                {
                    BesseljSol[i][n]=Besselj(n,k*Rarray[i]);
                    BesselySol[i][n]=Bessely(n,k*Rarray[i]);
                }
            }

            omp_set_num_threads(threads);
            #pragma omp parallel for shared(lmax, LegendreSol) //Generate Legendre Functions Recursively in Parallel for 0 to 180 degrees in integer increments
            for(int i=0;i<=180;i++)
            {
                for (int n=0; n<=lmax; ++n)
                {
                    LegendreSol[i][n]=LegendreP(n,cos( (double)i*PI/180.0 ));
                }
            }

        std::vector<double > pp (181,0);
            pp=ForwardNumerov(Potential,LegendreSol,h,lmax,scale,Energy,k);


    //=================================================================================================================================================================
    //    This is the end of the driving function and beginning of writing out data to file and/or console
    //=================================================================================================================================================================

    std::ofstream Solutions_out ("Scattering.txt");
      if (Solutions_out.is_open())
      {
            for (int b=0; b<=180; ++b)
            {
                //for (int c=0; c<=lmax; ++c)
                //{
                    Solutions_out<<pp[b]<<"\t";
                //}
                Solutions_out<<"\n";
            }
            Solutions_out.close();
      }

    //Timekeeping------------------------------------------------------------------------------------------------------------------------------------------------------
        endTOT=std::chrono::system_clock::now();
        std::chrono::duration<double> TOTtime=endTOT-startTOT;
        std::cout<<"\nTotal execution time: "<<TOTtime.count()<<std::endl<<std::endl<<"End of line, man."<<std::endl<<std::endl<<std::endl;
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------

    //Memory handling--------------------------------------------------------------------------------------------------------------------------------------------------
        std::vector<std::vector<double > >().swap(BesseljSol);
        std::vector<std::vector<double > >().swap(BesselySol);
        std::vector<std::vector<double > >().swap(LegendreSol);
        std::vector<std::vector<double > >().swap(Potential);
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------

    return 0;
}


//=================================================================================================================================================================
//    End of line, man
//=================================================================================================================================================================
