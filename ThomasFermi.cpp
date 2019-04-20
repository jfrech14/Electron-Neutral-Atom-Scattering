/* C++ Program by Josh Frechem
   Compile with g++ with OpenMP support or
   On MacOS:  clang++ -Xpreprocessor -fopenmp -lomp -o Scattering ThomasFermi.cpp    */

#include "Headers/Globals.h"

//Units are eV,A,s,
//Physical Constants---------------------------------------------------------------
    const double PI = 3.141592653589793238462643383279502884;
    const double h2m = 7.6359;                  //eV-A^2    
    const double SOL=299792458*pow(10,10);      //Speed of light in A/s
    const double m=0.5109989461*pow(10,6);      //electron mass eV/c^2    
    const double hbar=6.582119514*pow(10,-16);  //hbar eV-s
    const double hbar2m=h2m;//hbar*hbar*SOL*SOL/2/m;  //hbar^2/2m
    const double ee=14.409;                     //electron charge squared in eV-A
    const double aBohr=h2m/ee;                  //Bohr radius in A
    const double ch=1.6021766208*pow(10,-19);   //electron charge in C
//---------------------------------------------------------------------------------

//Setting up parallel processing-----------------------------------------------------------------------------------------------------------------------------------
    int numcores= omp_get_num_procs(); //set number of parallel workers by finding number of CPU workers. id1 is to be optimized for the amount of work per task
    double id1=1.0;
    const int threads=id1*(numcores);
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------


// Begin Simulation main function
int main()
{
    //Reading input file. Change the size of "inputs" if there are more inputs-------
    std::ifstream inputfile ("Input.txt");
    std::string line;
    std::vector<double> inputs (4,0.0);
    if (inputfile.is_open())
    {
        int inputi=0;
        while (std::getline(inputfile, line))
        {
            if (line != "#" )
            {
                std::istringstream iss(line);
                inputs[inputi/2]; 
                while ((iss >> inputs[inputi/2])){}
                inputi++;
            }
        }
    }

    else
    {
        std::cout << "Unable to open file";
        return 0;
    }
    //-------------------------------------------------------------------------------

    int MeshSize=inputs[3]; //mesh size
    int Z=inputs[0];           //atomic number
    double Energy=inputs[1];    //max energy level in eV to calculate
    double Rmax=2.0; //Mesh radius limit in angstroms
    double scale=5.0;            //Scale for how far out past Rmax to go
    double h=scale*Rmax/MeshSize;    //step size in angstroms
    double k=sqrt(Energy/hbar2m);                 //wavenumber 1/A
    double cval=-(2*Rmax)*(2*Rmax)*Energy/hbar2m; //constant for quadratic solving
    int lmax=(-1+sqrt(1-4*cval))/(2); //Set max lmax for highest energy for energy for loop
    int NumEnergies=inputs[2];
    
    //Timekeeping------------------------------------------------------------------------------------------------------------------------------------------------------
        std::chrono::time_point<std::chrono::system_clock> startTTM, endTTM, startDIFF,endDIFF, startTOT, endTOT;
        startTOT=std::chrono::system_clock::now();
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------



        std::cout<<"\n\n\n             Lenz-Jensen Atomic Potential Scattering\n   Shared Memory Parallelism Enabled:\t"<<numcores<<" processing cores used"<<std::endl;
        std::cout<<"__________________________________________________________________\n"<<std::endl;
        std::cout<<"Electron scattering off ";
        periodicsearch(Z);
        std::cout<<"\nNumber of partial waves to use: "<<lmax<<std::endl<<std::endl;

    //Memory Allocation------------------------------------------------------------------------------------------------------------------------------------------------
        std::vector<std::vector<double > > Potential (MeshSize+1,std::vector<double>(2,0));
        std::vector<double > Rarray (MeshSize,0.0);
        std::vector<double > NumerovSol (181,0);
        std::vector<std::vector<double > > DiffCrossSections (181,std::vector<double>(NumEnergies+1,0.0));
        std::vector< std::vector<double > > CrossSections (NumEnergies, std::vector<double> (2,0.0));
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------


    //=================================================================================================================================================================
    //   This is the end of simulation setup and beginning of driving function
    //=================================================================================================================================================================
    omp_set_num_threads(threads);
    #pragma omp parallel for shared(MeshSize, Rmax, Potential, Rarray,Z,scale)     
    for (int i=0;i<=MeshSize;++i)
    {
        Rarray[i]=scale*Rmax/MeshSize*i;  //First column is argument array
        Potential[i][0]=scale*Rmax/MeshSize*i;   //First column is argument array
        Potential[i][1]=ThomasFermiPotential(Z,Rarray[i]); //Second column is V array
        //Potential[i][1]=HarmonicPotential(Potential[i][0],1.0,1.0);
        //Potential[i][1]=ZeroPotential(Z,Potential[i][0]);
    }

    for (int o=0;o<=180;o++)
    {
        DiffCrossSections[o][0]=o;
    }

    omp_set_num_threads(threads);
    #pragma omp parallel for shared(CrossSections,lmax,scale,h,Potential,NumEnergies)
    for (int en=0;en<NumEnergies;en++)
    {
        CrossSections[en][0]=(double) ((en+1.0)*1000.0/NumEnergies);
        Energy=CrossSections[en][0];
        k=sqrt(Energy/hbar2m);  


        NumerovSol=ForwardNumerov(CrossSections,Potential,h,lmax,scale,Energy,k,en);
        for (int jj=0;jj<=180;jj++)
        {
            DiffCrossSections[jj][en+1]=NumerovSol[jj];
        }
    }
    //=================================================================================================================================================================
    //    This is the end of the driving function and beginning of writing out data to file and/or console
    //=================================================================================================================================================================

    std::ofstream Solutions1_out ("TotalCS.txt");
    if (Solutions1_out.is_open())
    {
        for (int b=0; b<NumEnergies; ++b)
        {
            Solutions1_out<<CrossSections[b][0];
            for (int c=1; c<=1; ++c)
            {
                Solutions1_out<<"\t"<<CrossSections[b][c];
                //std::cout<<"Total Cross Section for "<<CrossSections[b][0]<<"eV electron = "<<CrossSections[b][c]<<"  A^2\n"<<std::endl;
            }Solutions1_out<<std::endl;
        }
        Solutions1_out.close();
    }

    std::ofstream Solutions2_out ("DifferentialCS.txt");
    if (Solutions2_out.is_open())
    {
        for (int b=0; b<=180; ++b)
        {
            Solutions2_out<<DiffCrossSections[b][0];
            for (int c=1; c<=NumEnergies; ++c)
            {
                Solutions2_out<<"\t"<<DiffCrossSections[b][c];
            }Solutions2_out<<std::endl;
        }
        Solutions2_out.close();
    }

    //Timekeeping------------------------------------------------------------------------------------------------------------------------------------------------------
        endTOT=std::chrono::system_clock::now();
        std::chrono::duration<double> TOTtime=endTOT-startTOT;
        std::cout<<"\nTotal execution time: "<<TOTtime.count()<<std::endl<<std::endl<<"End of line, man."<<std::endl<<std::endl<<std::endl;
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------

    return 0;
}

//=================================================================================================================================================================
//    End of line, man
//=================================================================================================================================================================