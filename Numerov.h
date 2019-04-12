#ifndef Numerov_H
#define Numerov_H

#include "Globals.h"
#include "Potentials.h"


	double g(double V, double E)
	{
	    return(2/h2m*(V-E));
	}

	double numerov(int numE, int Mesh, int dh, const std::vector<std::vector<double > > &Potential)
	{
		std::vector<double> x(Mesh+1,0.0);
		std::vector<double> f(Mesh+1,0.0);
		std::vector<double> y(Mesh+1,0.0);
		std::vector<double> ysquared(Mesh+1,0.0);
		std::vector<double> V(Mesh+1,0.0);

		for(int i=0;i<=Mesh;++i)
		{
			//V[i]=Potential[i][1];
		}


		int s,iters,icl,zzero;
		double norm,Eupper,Elower,yicl,discont;


		int zero=numE;
		int nzeros=0;
		double Ecalc=0;


		//Begin loop to find Eigenvalue and integrate to find Eigenfunction
		for(;;)
		{
		    //Find Energy Eigenvalue and outward integration
		    Eupper=V[Mesh];
		    Elower=Eupper;
		    for(int i=0;i<=Mesh;i++)
		    {
		        if(V[i]<Elower)
		        {
		            Elower=V[i];
		        }
		        else{}
		        if(V[i]>Eupper)
		        {
		            Eupper=V[i];
		        }
		        else{}
		    }

		    if(Ecalc==0)
		    {
		        Ecalc=0.5*(Eupper+Elower);
		        iters=500;
		    }
		    else
		    {
		        iters=1;
		        std::cout<<"Here"<<std::endl;
		    }

		    for(int k=0;(k<iters) && (Eupper-Elower>pow(10,-10));k++)
		    {
		        f[0]=(dh*dh)/(12.)*g(V[0],Ecalc);
		        icl=-1;
		        for(int d=1;d<=Mesh;d++)
		        {
		            f[d]=(dh*dh)/(12.)*g(V[d],Ecalc);
		            if(f[d]==0)
		            {
		                f[d]=pow(10,-10);
		            }
		            else{}
		            if(f[d]!=std::copysign(f[d],f[d-1]))
		            {
		                icl=d;
		            }
		            else{}
		        }



		        for(int d=0;d<=Mesh;d++)
		        {
		            f[d]=1.0-f[d];
		            y[d]=0.0;
		        }

		        zzero=zero/2;

		        if(2*zzero==zero)
		        {
		            y[0]=1.0;
		            y[1]=0.5*(12.0-f[0]*10.0)*y[0]/f[1];
		        }
		        else
		        {
		            y[0]=0.0;
		            y[1]=1.0;
		        }

		        nzeros=0;
		        for(int a=1;a<=icl-1;a++)
		        {
		            y[a+1]=((12.-f[a]*10.)*y[a]-f[a-1]*y[a-1])/f[a+1];
		            if(y[a]!=std::copysign(y[a],y[a+1]))
		            {
		                ++nzeros;
		            }
		            else{}
		        }
		        yicl=y[icl];


		        if(2*zzero==zero)
		        {
		            nzeros=2*nzeros;
		        }
		        else
		        {
		            nzeros=2*nzeros+1;
		        }

		        if(iters>1)
		        {
		            if(nzeros!=zero)
		            {
		                if(k==0)
		                {
		                    std::cout<<"Iteration"<<"\t"<<"Energy"<<"\t\t"<<"Zeros\t\t"<<"Discont"<<"\n";
		                }
		                else{}

		                std::cout<<k<<"\t\t"<<std::setprecision(5)<<Ecalc<<"\t\t"<<nzeros<<"\t\t"<<discont<<"\n";

		                if(nzeros>zero)
		                {
		                    Eupper=Ecalc;
		                }
		                else
		                {
		                    Elower=Ecalc;
		                }
		                Ecalc=0.5*(Eupper+Elower);
		            }
		            else{}
		        }
		        else{}

		        //Begin inward integration after finding correct number of zeros
		        if(iters==1 || nzeros==zero)
		        {
		            y[Mesh]=dh;
		            y[Mesh-1]=(12.-10.*f[Mesh])*y[Mesh]/f[Mesh-1];
		            y[Mesh+1]=0;
		            for(int b=Mesh-1;b>=icl+1;b--)
		            {
		                y[b-1]=((12.-10.*f[b])*y[b]-f[b+1]*y[b+1])/f[b-1];
		            }

		            yicl /= y[icl];
		            for(int b=icl;b<=Mesh;b++)
		            {
		                y[b]*=yicl;
		            }

		            norm=0;
		            for(int b=1;b<=Mesh;b++)
		            {
		                norm+=(y[b]*y[b]);
		            }
		            norm=dh*(2.0*norm+(y[0]*y[0]));
		            norm=-sqrt(norm);

		            for(int b=0;b<=Mesh;b++)
		            {
		                y[b]/=norm;
		            }

		            if(iters>1)
		            {
		                s=icl;
		                discont=(y[s+1]+y[s-1]-(14.-12.*f[s])*y[s])/dh;
		                std::cout<<k<<"\t\t"<<Ecalc<<"\t\t"<<nzeros<<"\t\t"<<discont<<"\n";
		                if(discont*y[s]>0.0)
		                {
		                    Eupper=Ecalc;
		                }
		                else
		                {
		                    Elower=Ecalc;
		                }
		                Ecalc=0.5*(Eupper+Elower);
		            }

		        }

		    }
		if(nzeros==zero)
		        break;
		}
		return (Ecalc);
	}
#endif
