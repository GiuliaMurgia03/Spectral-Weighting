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
    fitsfile *fptr;

    // Open fits file
    if (!sp.open(argv[1], &fptr))
    {
        return EXIT_FAILURE;
    }

    // Splat
    // if (!sp.splat(fptr, "splat.fits"))
    //  {
    //      return EXIT_FAILURE;
    //  }

    // Smoothing
    if (!sp.gaussian_smoothing(fptr, "smooth.fits", 3))
    {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}