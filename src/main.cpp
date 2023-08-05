#include <iostream>
#include "SpectralWeighting.h"

using namespace std;
using namespace spacew;

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cout << "Usage: SpectralWeighting input.txt" << endl;
        return EXIT_FAILURE;
    }

    SpectralWeighting sp;

    // Splat
    if (!sp.splat(argv[1], "splat.fits"))
    {
        return EXIT_FAILURE;
    }

    // Smoothing
    // if (!sp.gaussian_smoothing(fptr, "smooth.fits", 3))
    //{
    //    return EXIT_FAILURE;
    //}

    return EXIT_SUCCESS;
}