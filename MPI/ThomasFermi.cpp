/* C++ Program by Josh Frechem
   Compile with g++ with MPI and OpenMP support or
   On MacOS:  clang++ -Xpreprocessor -fopenmp -lomp -o Scattering ThomasFermi.cpp
   Make the number of energies a multiple of how many processors to use*/

#include "Headers/Globals.h"
#define  MASTER		0

//Units are eV,A,s,
//Physical Constants---------------------------------------------------------------
    const double PI = 3.141592653589793238462643383279502884;
    const double h2m = 7.6359;                  //eV-A^2
    const double SOL=299792458*pow(10,10);      //Speed of light in A/s
    const double m=0.5109989461*pow(10,6);      //electron mass eV/c^2
    const double hbar=6.582119514*pow(10,-16);  //hbar eV-s
    const double hbar2m=hbar*hbar*SOL*SOL/2.0/m;//hbar^2/2m
    const double ee=14.409;                     //electron charge squared in eV-A
    const double aBohr=h2m/ee;                  //Bohr radius in A
    const double ch=1.6021766208*pow(10,-19);   //electron charge in C
//---------------------------------------------------------------------------------


// Begin Simulation main function
int main()
{
    int tag1, tag2, tag3, tag4,tag5,tag6, dest,offset,source,chunksize;
    //Reading input file. Change the size of "inputs" if there are more inputs-------
    std::ifstream inputfile ("Input.txt");
    std::string line;
    std::vector<double> inputs (6,0.0);
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

    int MeshSize=inputs[3];     //mesh size
    int Z=inputs[0];            //atomic number
    double Energy=inputs[1];    //max energy level in eV to calculate
    double Rmax=2.0;            //Mesh radius limit in angstroms
    double scale=inputs[4];           //Scale for how far out past Rmax to go
    double h=scale*Rmax/MeshSize;                 //step size in angstroms
    double k=sqrt(Energy/hbar2m);                 //wavenumber 1/A
    double cval=-(2*Rmax)*(2*Rmax)*Energy/hbar2m; //constant for quadratic solving
    int lmax=(-1+sqrt(1-4*cval))/(2);             //Set max lmax for highest energy for energy for loop
    int NumEnergies=inputs[2]; //Number of energies to calculate the cross section for
    double Espec=inputs[5];    //Specific energy to get partial waves at

    //Memory Allocation------------------------------------------------------------------------------------------------------------------------------------------------
        std::vector<std::vector<double > > Potential (MeshSize+1,std::vector<double>(2,0));
        std::vector<double > Rarray (MeshSize+1,0.0);
        std::vector<double > NumerovSol (181,0);
        std::vector<std::vector<double > > DiffCrossSections (NumEnergies+1,std::vector<double>(181,0.0));
        std::vector< std::vector<double > > CrossSections (2, std::vector<double> (NumEnergies,0.0));
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------


    //=================================================================================================================================================================
    //   This is the end of simulation setup and beginning of driving function
    //=================================================================================================================================================================

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
        DiffCrossSections[o][0]=(double) o;
    }

    for (int en=0;en<NumEnergies;en++)
    {
        CrossSections[0][en]=(double) pow(10,(-2.0+(en)*(2+log10(inputs[1]))/NumEnergies));
    }

    MPI_Status status;
    MPI_Init(NULL,NULL);
    int rank,size,name_len;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name, &name_len);
    tag2=1;
    tag1=2;
    tag3=4;
    chunksize = (NumEnergies / size);


    if (rank == MASTER)
    {
        //Master Timekeeping------------------------------------------------------------------------------------------------------------------------------------------------------
            std::chrono::time_point<std::chrono::system_clock> startTTM, endTTM, startDIFF,endDIFF, startTOT, endTOT;
            startTOT=std::chrono::system_clock::now();
        //-----------------------------------------------------------------------------------------------------------------------------------------------------------------

            std::cout<<"\n\n\n             Lenz-Jensen Atomic Potential Scattering\n   MPI Parallel Computing Enabled:\t"<<size<<" processing cores used"<<std::endl;
            std::cout<<"__________________________________________________________________\n"<<std::endl;
            std::cout<<"Electron scattering off "<<periodicsearch(Z);
            std::cout<<"\nNumber of partial waves to use: "<<lmax<<std::endl<<std::endl;

        /* Send each task its portion of the array - master keeps 1st part */
        offset = chunksize;
        for (dest=1; dest<size; dest++)
        {
            MPI_Send(&offset, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
            offset = offset + chunksize;
        }
        
        /* Master does its part of the work */
        offset = 0;
        std::cout<<"Processor: "<<processor_name<<"   rank:"<<rank<<"  Energies="<<CrossSections[0][offset]<<" - "<<CrossSections[0][offset+chunksize-1]<<std::endl;
        for(int i=offset; i < offset + chunksize; i++)
        {
            Energy=CrossSections[0][i];
            k=sqrt(Energy/hbar2m);
            NumerovSol=ForwardNumerov(Espec,CrossSections,Potential,h,lmax,scale,Energy,k,i);
            for (int jj=0;jj<=180;jj++)
            {
                DiffCrossSections[i+1][jj]=NumerovSol[jj];
            }
        }

        /* Wait to receive results from each slave task */
        for (int i=1; i<size; i++)
        {
            source = i;
            MPI_Recv(&offset, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
            MPI_Recv(&CrossSections[1][offset], chunksize, MPI_DOUBLE, source, tag2, MPI_COMM_WORLD, &status);
            MPI_Recv(&DiffCrossSections[offset+1][0], (181)*chunksize, MPI_DOUBLE, source, tag3, MPI_COMM_WORLD, &status);
        }

        //=================================================================================================================================================================
        //    This is the beginning of writing out data to file and/or console from Master
        //=================================================================================================================================================================
        std::ofstream Solutions1_out (periodicsearch(Z)+"TotalCS.txt");
        if (Solutions1_out.is_open())
        {
            for (int b=0; b<NumEnergies; ++b)
            {
                Solutions1_out<<CrossSections[0][b];
                for (int c=1; c<=1; ++c)
                {
                    Solutions1_out<<"\t"<<CrossSections[c][b];
                }Solutions1_out<<std::endl;
            }
            Solutions1_out.close();
        }

        std::ofstream Solutions2_out (periodicsearch(Z)+"DifferentialCS.txt");
        if (Solutions2_out.is_open())
        {
            for (int b=0; b<=180; ++b)
            {
                Solutions2_out<<DiffCrossSections[0][b];
                for (int c=1; c<NumEnergies+1; ++c)
                {
                    Solutions2_out<<"\t"<<DiffCrossSections[c][b];
                }Solutions2_out<<std::endl;
            }
            Solutions2_out.close();
        }

        std::ofstream Solutions3_out (periodicsearch(Z)+"Potential.txt");
        if (Solutions3_out.is_open())
        {
            for (int b=0; b<=MeshSize; ++b)
            {
                Solutions3_out<<Potential[b][0];
                for (int c=1; c<=1; ++c)
                {
                    Solutions3_out<<"\t"<<(1.0/hbar2m)*( Potential[b][c] + hbar2m*2.0/Rarray[b]/Rarray[b]);
                }Solutions3_out<<std::endl;
            }
            Solutions3_out.close();
        }

        //End of Master Timekeeping----------------------------------------------------------------------------------------------------------------------------------------
            endTOT=std::chrono::system_clock::now();
            std::chrono::duration<double> TOTtime=endTOT-startTOT;
            std::cout<<"\nTotal execution time: "<<TOTtime.count()<<std::endl<<std::endl<<"End of line, man."<<std::endl<<std::endl<<std::endl;
        //-----------------------------------------------------------------------------------------------------------------------------------------------------------------
    }


    //Non-master tasks
    if (rank > MASTER)
    {
        /* Receive my portion of array from the master task */
        source = MASTER;
        MPI_Recv(&offset, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);

        std::cout<<"Processor: "<<processor_name<<"   rank:"<<rank<<"  Energies="<<CrossSections[0][offset]<<" - "<<CrossSections[0][offset+chunksize-1]<<std::endl;

        /* Each task does its part of the work */
        for(int i=offset; i < offset + chunksize; i++)
        {
            Energy=CrossSections[0][i];
            k=sqrt(Energy/hbar2m);
            NumerovSol=ForwardNumerov(Espec,CrossSections,Potential,h,lmax,scale,Energy,k,i);
            for (int jj=0;jj<=180;jj++)
            {
                DiffCrossSections[i+1][jj]=NumerovSol[jj];
            }
        }

        /* Send task results back to the master task */
        dest = MASTER;
        MPI_Send(&offset, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
        MPI_Send(&CrossSections[1][offset], chunksize, MPI_DOUBLE, dest, tag2, MPI_COMM_WORLD);
        MPI_Send(&DiffCrossSections[offset+1][0], (181)*chunksize, MPI_DOUBLE, dest, tag3, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}

//=================================================================================================================================================================
//    End of line, man
//=================================================================================================================================================================