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
        bool open(const string &);
    };

}