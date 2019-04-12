#ifndef Potentials_H
#define Potentials_H

	double ZeroPotential(double Z, double r)
	{
		return(0);
	}
	double ThomasFermiPotential(double Z, double r)// V=eV, r=angstroms
	{
	    double rr=r*pow(10,10);
	    double b2=0.3344, b3=0.0485, b4=0.002647;
	    long double x=4.5397*pow(Z,1.0/6.0)*sqrt(rr);
	    long double V=-Z*ee/rr*exp(-1.0*x)*(1 + x + b2*x*x + b3*x*x*x + b4*x*x*x*x);
	    return(V);
	}

	double HarmonicPotential(double x, double omega, double m)
	{
	    return(m*omega*omega*0.5*x*x);
	}

#endif
