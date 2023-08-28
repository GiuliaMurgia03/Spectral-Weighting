#include <iostream>
#include "SpectralWeighting.h"
#include "SpectralSimulations.h"

using namespace std;
using namespace spacew;

int main(int argc, char *argv[])

{

    float sigma = 0.0;
    string option = "-help";
    if (argc >= 2)
    {
        option = argv[1];
    }

    // Check smooth option
    if (option == "-smooth" && argc >= 3)
    {
        sigma = std::stof(argv[2]);
        for (int i = 3; i < argc; i++)
        {
            argv[i - 2] = argv[i];
        }
        argc = argc - 2;
        option = argv[1];
    }

    cout << "Option :" << option << endl;

    // Help
    if (((argc == 2 && option == "-help") || argc <= 1))
    {
        cout << "Type one of the following:" << endl;
        cout << "SpectralWeighting -splat input.fits output.fits [bchan echan]" << endl;
        cout << "SpectralWeighting [-smooth sigma] -weighted_splat input.fits output.fits size [bchan echan]" << endl;
        cout << "SpectralWeighting -gaussian_smooth sigma input.fits output.fits" << endl;
        cout << "SpectralWeighting -local_noise input.fits output.fits size" << endl;
        cout << "SpectralWeighting [-smooth sigma] -local_weights input.fits output.fits size" << endl;
        cout << "SpectralWeighting -merge filelist.txt output.fits [bchan echan]" << endl;
        cout << "SpectralWeighting [-smooth sigma] -weighted_merge exponent filelist.txt output.fits size [bchan echan]" << endl;
        cout << "SpectralWeighting -simul noise_model outfile.fits nx ny nz nsources sigma_noise [rfi_infile]" << endl;
        cout << "SpectralWeighting -sum input1.fits input2.fits outfile.fits" << endl;
        cout << "Valid noise models are: white_noise, vertical_noise, horizontal_noise" << endl;

        return EXIT_SUCCESS;
    }

    SpectralWeighting sp;
    SpectralSimulations ssimul;

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
    if (argc == 5 && option == "-gaussian_smooth" && !sp.gaussian_smoothing(argv[3], argv[4], std::stof(argv[2])))
    {
        return EXIT_FAILURE;
    }

    // Sum
    if (argc == 5 && option == "-sum" && !sp.sum_fits(argv[2], argv[3], argv[4]))
    {
        return EXIT_FAILURE;
    }

    // Local Noise
    if (argc == 5 && option == "-local_noise" && !sp.local_noise(argv[2], argv[3], std::stoi(argv[4])))
    {
        return EXIT_FAILURE;
    }

    // Local Weight
    if (argc == 5 && option == "-local_weights" && !sp.local_weights(argv[2], argv[3], std::stoi(argv[4]), sigma))
    {
        return EXIT_FAILURE;
    }

    // Weighted Splat
    if (argc >= 4 && option == "-weighted_splat")
    {
        if (argc == 5 && !sp.weighted_splat(argv[2], argv[3], std::stoi(argv[4]), 0, 0, sigma))
        {
            return EXIT_FAILURE;
        }

        if (argc == 7 && !sp.weighted_splat(argv[2], argv[3], std::stoi(argv[4]), std::stoi(argv[5]), std::stoi(argv[6]), sigma))
        {
            return EXIT_FAILURE;
        }
    }

    // Weighted Merge
    if (argc >= 5 && option == "-weighted_merge")
    {
        if (argc == 6 && (!sp.set_exponent(std::stof(argv[2])) || !sp.weighted_merge(argv[3], argv[4], std::stoi(argv[5]), 0, 0, sigma)))
        {
            return EXIT_FAILURE;
        }

        if (argc == 8 && (!sp.set_exponent(std::stof(argv[2])) || !sp.weighted_merge(argv[3], argv[4], std::stoi(argv[5]), std::stoi(argv[6]), std::stoi(argv[7]), sigma)))

        {
            return EXIT_FAILURE;
        }
    }

    // Merge
    if (argc >= 3 && option == "-merge")
    {
        if (argc == 4 && !sp.weighted_merge(argv[2], argv[3], 0.0))
        {
            return EXIT_FAILURE;
        }

        if (argc == 6 && !sp.weighted_merge(argv[2], argv[3], 0.0, std::stoi(argv[4]), std::stoi(argv[5])))
        {
            return EXIT_FAILURE;
        }
    }

    // Simulation
    if (argc >= 9 && option == "-simul")
    {
        if (argc == 9)
        {
            if (!ssimul.set_noise_model(argv[2]) || !ssimul.create_outfits(argv[3], std::stoi(argv[4]), std::stoi(argv[5]), std::stoi(argv[6]), std::stoi(argv[7]), std::stof(argv[8])))

                return EXIT_FAILURE;
        }

        if (argc == 10)
        {
            if (!ssimul.set_noise_model(argv[2]) || !ssimul.create_outfits(argv[3], std::stoi(argv[4]), std::stoi(argv[5]), std::stoi(argv[6]), std::stoi(argv[7]), std::stof(argv[8]), argv[9]))

                return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}