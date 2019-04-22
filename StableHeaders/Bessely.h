#ifndef Bessely_H
#define Bessely_H

    double Bessely(int l, double x)
    {
        double y2=0.0;
        double y0=-cos(x)/x;
        double y1=y0/x-sin(x)/x;
        if(l==0)
            return y0;
        if(l==1)
            return y1;

        for(int i=2;i<=l;++i)
        {
            y2=(2.0*i-1)/x*y1-y0;
            y0=y1;
            y1=y2;
        }
        return(y2);
    }

#endif
