#ifndef LegendreP_H
#define LegendreP_H

	double LegendreP(int l, double x)
	{
        double j0=1;
        if (l<=0) return j0;
        double j1=x;
        double j2=j1;
        for (int i=2; i<=l; i++)
        {
            j2=j1*(2*i-1)*x/(i) - (i-1)/(i)*j0;
            j0=j1;
            j1=j2;
        }
        return j2;
	}

#endif