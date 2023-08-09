#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include "SpectralSimulations.h"
#include "fits.h"

using namespace std;

namespace spacew
{

    SpectralSimulations::SpectralSimulations()
    {
    }

    SpectralSimulations::~SpectralSimulations()
    {
    }

    bool SpectralSimulations::init_flux_bins(const float &s_min, const float &s_max, const float &delta_s)
    {

        min_flux = s_min;
        max_flux = s_max;
        ds = delta_s;
        int nbins = (max_flux - min_flux) / ds;
        flux_bins.resize(nbins);

        // Fill flux vector according to logN-logS (dN/dS=c*S^-5/2)

        flux_bins[0] = pow(0.5 * ds, -5.0 / 2.0) * ds;

        for (int i = 1; i < nbins; i++)
        {
            flux_bins[i] = flux_bins[i - 1] + pow((0.5 + i) * ds, -5.0 / 2.0) * ds;
        }

        // Normalize dividing by the last flux bin so that the cumulative distribution tends to 1
        for (int i = 0; i < nbins; i++)
        {
            flux_bins[i] = flux_bins[i] / flux_bins.back();
        }
    }

    bool SpectralSimulations::create_outfits(const string &outfile, int nx, int ny, int nz, int nsources, float sigma_noise)
    {

        spacew::fits outfits;
        if (!outfits.create(outfile))
        {
            return false;
        }

        int status = 0;
        int bitpix = -32;
        int naxis = 3;
        long naxes[naxis];
        naxes[0] = nx;
        naxes[1] = ny;
        naxes[2] = nz;
        fits_create_img(outfits.get_fptr(), bitpix, naxis, naxes, &status);
        outfits.get_header_info();

        if (status) // Check that worked
        {
            fits_report_error(stderr, status);
            return false;
        }

        outfits.fill(0);
        add_pointsources_model(outfits, nsources,sigma_noise);

        if (!outfits.close()) // Check that worked
        {
            return false;
        }

        return true;
    }

    bool SpectralSimulations::add_pointsources_model(fits &model_fits, int nsources, float sigma_noise)
    {
        init_flux_bins(0.001, 1, 0.00001);

        int nx = model_fits.get_naxes(0);
        int ny = model_fits.get_naxes(1);
        int nz = model_fits.get_naxes(2);

        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_int_distribution<uint32_t> uniform_random_nx(1, nx);
        std::uniform_int_distribution<uint32_t> uniform_random_ny(1, ny);

        std::normal_distribution<float> gauss_noise(0, sigma_noise);

        rng.seed(1); // Always the same random simulation

        vector<int> xpos(nsources);
        vector<int> ypos(nsources);
        vector<float> image(nx * ny);

        for (int n = 0; n < nsources; n++)
        {
            xpos[n] = uniform_random_nx(rng);
            ypos[n] = uniform_random_ny(rng);
        }

        for (int k = 0; k < nz; k++)
        {
            cout << "Working on channel: " << k + 1 << " of " << nz << "\t\r" << std::flush;
            model_fits.read_channel_image(k, image);

            for (int n = 0; n < nsources; n++)
            {

                image[xpos[n] + nx * ypos[n]] += extract_random_flux();
            }

            for (int i = 0; i < nx * ny; i++)
            {
                image[i] += gauss_noise(rng);
            }

            model_fits.write_channel_image(k, image);
        }

        return true;
    }

    float SpectralSimulations::extract_random_flux()
    {

        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_real_distribution<> random_number(0.0, 1.0);

        float rn = random_number(rng);

        for (int i = 0; i < flux_bins.size(); i++)
        {
            // cout<<i<<" "<<flux_bins[i]<<" "<<rn<<endl;

            if (flux_bins[i] > rn)
            {
                return min_flux + (i + 0.5) * ds;
            }
        }

        return 1;
    }

}
