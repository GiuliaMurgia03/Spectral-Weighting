#pragma once
#include <string>

using namespace std;

namespace spacew
{

    class rfi
    {
    private:
        int xpix1, xpix2, ypix1, ypix2;
        float intensity;
        int bchan, echan;    
        float m, q;
    public:
        rfi(int x1, int y1, int x2, int y2, float i, int ch1, int ch2);
        ~rfi();
        float get_intensity(int xp, int yp, int ch);

    };

}