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
        bool splat(const string& infile, const string& outfile);
        bool gaussian_smoothing(fitsfile*, const string&, float sigma);
    };

}