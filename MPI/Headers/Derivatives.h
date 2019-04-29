#ifndef Derivatives_H
#define Derivatives_H

//Takes an array and a step size. This method was chosen as more useful than a 2D vector with X and Y data
    
    std::vector<double> FirstDerivative(std::vector<double>& Variable, double dh)  
    {
        int size=Variable.size();
        std::vector<double> d1(size,0.0);
        for (int i=0;i<=size;++i)
        {
            if (i<2)  //forward finite difference
            {
                d1[i]=(Variable[i+1]-Variable[i])/dh;
            }
            else if (i>=2&&i<=size-3)  //central finite difference
            {
                d1[i]=(1.0/12.0/dh)*(Variable[i-2]-8.0*Variable[i-1]+8.0*Variable[i+1]-Variable[i+2]);
            }
            else //backward finite difference
            {
                d1[i]=(Variable[i]-Variable[i-1])/dh;
            }
        }
        return(d1);
    }

    std::vector<double> SecondDerivative(std::vector<double>& Variable, double dh)  
    {
        int size=Variable.size();
        std::vector<double> d1(size,0.0);
        for (int i=0;i<=size;++i)
        {
            if (i<2)  //forward finite difference
            {
                d1[i]=(Variable[i+2] -2.0*Variable[i+1]+Variable[i])/dh/dh;
            }
            else if (i>=2&&i<=size-3)  //central finite difference
            {
                d1[i]=(Variable[i-1]-2.0*Variable[i]+Variable[i+1])/dh/dh;
            }
            else //backward finite difference
            {
                d1[i]=(Variable[i]-2.0*Variable[i-1]+Variable[i-2])/dh/dh;
            }
        }
        return(d1);
    }

#endif