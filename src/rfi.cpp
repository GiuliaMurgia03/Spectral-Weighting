#include <iostream>
#include "rfi.h"


using namespace std;

namespace spacew
{

    rfi::rfi(int x1, int y1, int x2, int y2, float i, int ch1, int ch2)
    {
        xpix1=x1;
        xpix2=x2;
        ypix1=y1;
        ypix2=y2;
        intensity=i;
        bchan=ch1; 
        echan=ch2;   
        if(x1!=x2){
            m=(1.0*y2-y1)/(x2-x1);
            q=y1-(m*x1);
        }
        else
        {
            m=0;
            q=0;
        }
        
    }

    rfi::~rfi()
    {
    }

    float rfi::get_intensity(int xp, int yp, int ch) {
        if(ch<bchan || ch>echan){
            return 0.0;
        }


    if(m !=0 && (xp<std::min(xpix1, xpix2) || xp > std::max(xpix1, xpix2)))
    {
        return 0.0; 
    }
    if(m ==0 && (yp<std::min(ypix1, ypix2) || yp > std::max(ypix1, ypix2)))
    {
        return 0.0; 
    }
    

        float d=0;

        if(m==0 && q==0) {
            d=fabs(xp-xpix1);
        }
        else{
            d=(fabs(yp-(m*xp+q)))/sqrt(1+pow(m,2));
        }

        if(d<1.5){
            return intensity; 
        }
    

        return 0.0;
    }

}