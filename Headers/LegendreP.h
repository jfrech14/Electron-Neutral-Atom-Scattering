#ifndef LegendreP_H
#define LegendreP_H

	double LegendreP(int l, double x)
	{
        double P0=1.0;
        if (l<=0) 
            {
                return P0;
            }
        double P1=x;
        if (l==1) 
            {
                return P1;
            }
        double P2=P1;
        for (int i=2; i<=l; ++i)
        {
            P2=P1/i*(2.0*i-1)*x - (double)(i-1)/(i)*P0;
            P0=P1;
            P1=P2;
        }
        return P2;
	}

#endif