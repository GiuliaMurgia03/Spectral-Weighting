#include <iostream>
#include "SpectralWeighting.h"

using namespace std;
using namespace spacew;

int main(int argc, char *argv[])
{

    string option = argv[1];
    cout << "Option :" << option << endl;

    // Help
    if ((argc == 2 && option == "-help" )|| argc != 4)
    {
        cout << "Type one of the following:" << endl;
        cout << "SpectralWeighting -splat input.fits output.fits" << endl;
        cout << "SpectralWeighting -smooth input.fits output.fits" << endl;
        return EXIT_SUCCESS;
    }

    SpectralWeighting sp;

    // Splat
    if (argc == 4 && option == "-splat" && !sp.splat(argv[2], argv[3]))
    {
        return EXIT_FAILURE;
    }

    // Smooth
    if (argc == 4 && option == "-smooth" && !sp.gaussian_smoothing(argv[2], argv[3], 3))
    {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}