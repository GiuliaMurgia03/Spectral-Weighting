#pragma once
#include "fitsio.h"
#include <string>

using namespace std;

namespace spacew
{

    class SpectralWeighting
    {
    public:
        SpectralWeighting();
        ~SpectralWeighting();
        bool test(const string&);
        bool open(const string&, fitsfile**);
        bool splat(fitsfile*, const string&);
    };

}