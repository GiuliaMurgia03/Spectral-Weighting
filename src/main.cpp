#include <iostream>
#include "SpectralWeighting.h"

using namespace std;
using namespace spacew;

int main(int argc, char *argv[])
{

    string option = argv[1];
    cout << "Option :" << option << endl;

    // Help
    if ((argc == 2 && option == "-help"))
    {
        cout << "Type one of the following:" << endl;
        cout << "SpectralWeighting -splat input.fits output.fits [bchan echan]" << endl;
        cout << "SpectralWeighting -weighted_splat input.fits output.fits size [bchan echan]" << endl;
        cout << "SpectralWeighting -smooth input.fits output.fits" << endl;
        cout << "SpectralWeighting -local_noise input.fits output.fits size" << endl;
        cout << "SpectralWeighting -local_weights input.fits output.fits size" << endl;
        cout << "SpectralWeighting -weighted_merge filelist.txt output.fits size [bchan echan]" << endl;

        return EXIT_SUCCESS;
    }

    SpectralWeighting sp;

    // Splat
    if (argc >= 4 && option == "-splat")
    {
        if (argc == 4 && !sp.splat(argv[2], argv[3]))
        {
            return EXIT_FAILURE;
        }

        if (argc == 6 && !sp.splat(argv[2], argv[3], std::stoi(argv[4]), std::stoi(argv[5])))
        {
            return EXIT_FAILURE;
        }
    }

    // Smooth
    if (argc == 4 && option == "-smooth" && !sp.gaussian_smoothing(argv[2], argv[3], 3))
    {
        return EXIT_FAILURE;
    }

    // Local Noise
    if (argc == 5 && option == "-local_noise" && !sp.local_noise(argv[2], argv[3], std::stoi(argv[4])))
    {
        return EXIT_FAILURE;
    }

    // Local Weight
    if (argc == 5 && option == "-local_weights" && !sp.local_weights(argv[2], argv[3], std::stoi(argv[4])))
    {
        return EXIT_FAILURE;
    }

    // Weighted Splat
    if (argc >= 4 && option == "-weighted_splat")
    {
        if (argc == 5 && !sp.weighted_splat(argv[2], argv[3], std::stoi(argv[4])))
        {
            return EXIT_FAILURE;
        }

        if (argc == 7 && !sp.weighted_splat(argv[2], argv[3], std::stoi(argv[4]), std::stoi(argv[5]), std::stoi(argv[6])))
        {
            return EXIT_FAILURE;
        }
    }


    // Weighted Merge
    if (argc >= 4 && option == "-weighted_merge")
    {
        if (argc == 5 && !sp.weighted_merge(argv[2], argv[3], std::stoi(argv[4])))
        {
            return EXIT_FAILURE;
        }

        if (argc == 7 && !sp.weighted_merge(argv[2], argv[3], std::stoi(argv[4]), std::stoi(argv[5]), std::stoi(argv[6])))
        {
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}