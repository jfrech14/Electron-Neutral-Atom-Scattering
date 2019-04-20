#ifndef Potentials_H
#define Potentials_H

	double ZeroPotential(double Z, double r)
	{
		return(0);
	}

	double ThomasFermiPotential(double Z, double rr)// V=eV, r=A
	{
	    //double rr=r*pow(10,10);//converting to angstroms for this particular equation
	    double b2=0.3344, b3=0.0485, b4=0.002647;
	    double x=4.5397*pow(Z,1.0/6.0)*sqrt(rr);
	    double V=-Z*ee/rr*exp(-1.0*x)*(1.0 + x + b2*x*x + b3*x*x*x + b4*x*x*x*x);
	    return(V);
	}

	double HarmonicPotential(double x, double omega, double m)
	{
	    return(m*omega*omega*0.5*x*x);
	}

	double SquareWellPotential(double x, double rwell)
	{
	    double FSWpotential=-5.0; //Attractive potential if negative
	    if (x<=rwell)
	    {
	    	return(FSWpotential);
	    }
	    else
	    {
	    	return (0);
	    }
	}

#endif
