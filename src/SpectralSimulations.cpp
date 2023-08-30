#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <sstream>
#include <chrono>
#include "SpectralSimulations.h"
#include "fits.h"

using namespace std;

namespace spacew
{

    SpectralSimulations::SpectralSimulations()
    {
        noise_model = "white_noise";
    }

    SpectralSimulations::~SpectralSimulations()
    {
    }

    bool SpectralSimulations::set_noise_model(const string &model)
    {
        if (model != "white_noise" && model != "vertical_noise" && model != "horizontal_noise")
        {
            cout << "Model noise not implemented" << endl;
            return false;
        }
        noise_model = model;
        return true;
    }

    bool SpectralSimulations::init_flux_bins(const float &s_min, const float &s_max, const float &delta_s)
    {
        // Sample in log scale

        min_flux = s_min;
        max_flux = s_max;
        ds = delta_s;
        int nbins = (max_flux - min_flux) / ds;
        flux_bins.resize(nbins);

        // Fill flux vector according to logN-logS (dN/dS=c*S^-5/2)

        float flux = min_flux + (0.5) * ds;

        flux_bins[0] = pow(flux, -5.0 / 2.0) * ds;

        for (int i = 1; i < nbins; i++)
        {
            flux = min_flux + (0.5 + i) * ds;
            flux_bins[i] = flux_bins[i - 1] + pow(flux, -5.0 / 2.0) * ds;
        }

        // Normalize dividing by the last flux bin so that the cumulative distribution tends to 1
        for (int i = 0; i < nbins; i++)
        {
            flux_bins[i] = flux_bins[i] / flux_bins.back();
        }

        return true;
    }

    bool SpectralSimulations::create_outfits(const string &outfile, int nx, int ny, int nz, int nsources, float sigma_noise, string rfi_infile)
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

        if (nsources > 0)
        {
            add_pointsources_model(outfits, nsources);
        }
        if (sigma_noise > 0)
        {
            if (noise_model == "white_noise")
            {
                add_noise_model(outfits, sigma_noise);
            }
            else if (noise_model == "vertical_noise")
            {
                add_vertical_scan_noise_model(outfits, sigma_noise, 0.8);
            }

            else if (noise_model == "horizontal_noise")
            {
                add_horizontal_scan_noise_model(outfits, sigma_noise, 0.8);
            }
        }

        if (rfi_infile != "")
        {
            add_RFI(outfits, rfi_infile);
        }

        if (!outfits.close()) // Check that worked
        {
            return false;
        }

        return true;
    }

    bool SpectralSimulations::add_pointsources_model(fits &model_fits, int nsources)
    {
        init_flux_bins(1e-4, 1, 1e-4);

        int nx = model_fits.get_naxes(0);
        int ny = model_fits.get_naxes(1);
        int nz = model_fits.get_naxes(2);

        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_int_distribution<uint32_t> uniform_random_nx(1, nx);
        std::uniform_int_distribution<uint32_t> uniform_random_ny(1, ny);

        rng.seed(1); // Always the same random simulation

        vector<int> xpos(nsources);
        vector<int> ypos(nsources);
        vector<float> image(nx * ny);
        vector<float> sources_image(nx * ny);

        for (int n = 0; n < nsources; n++)
        {
            xpos[n] = uniform_random_nx(rng);
            ypos[n] = uniform_random_ny(rng);
            sources_image[xpos[n] + nx * ypos[n]] += extract_random_flux();
        }

        for (int k = 0; k < nz; k++)
        {
            cout << "Working on channel: " << k + 1 << " of " << nz << "\t\r" << std::flush;
            model_fits.read_channel_image(k, image);

            for (int idx = 0; idx <= nx * ny; idx++)
            {
                image[idx] += sources_image[idx];
            }

            model_fits.write_channel_image(k, image);
        }

        return true;
    }

    bool SpectralSimulations::add_spectral_line_model(fits &model_fits, int xp, int yp, int channel, float peak_intensity, float width)
    {
        int nx = model_fits.get_naxes(0);
        int ny = model_fits.get_naxes(1);
        int nz = model_fits.get_naxes(2);

        vector<float> image(nx * ny);

        for (int k = 0; k < nz; k++)
        {
            cout << "Working on channel: " << k + 1 << " of " << nz << "\t\r" << std::flush;
            model_fits.read_channel_image(k, image);

            image[xp + yp * nx] += peak_intensity * exp(-0.5 * pow((channel - k) / width, 2.0));

            model_fits.write_channel_image(k, image);
        }

        return true;
    }

    bool SpectralSimulations::add_noise_model(fits &model_fits, float sigma_noise)
    {
        // Get a random seed from the clock
        // Example from https://cplusplus.com/reference/random/mersenne_twister_engine/seed/
        typedef std::chrono::high_resolution_clock myclock;
        myclock::time_point beginning = myclock::now();

        int nx = model_fits.get_naxes(0);
        int ny = model_fits.get_naxes(1);
        int nz = model_fits.get_naxes(2);

        std::random_device dev;
        std::mt19937 rng(dev());
        std::normal_distribution<float> gauss_noise(0, sigma_noise);

        // Always the same random simulation
        // rng.seed(1);

        // Obtain a seed from the timer
        myclock::duration d = myclock::now() - beginning;
        unsigned s = d.count();
        rng.seed(s);
        cout << "Using random seed: " << s << endl;

        vector<float> image(nx * ny);

        for (int k = 0; k < nz; k++)
        {
            cout << "Working on channel: " << k + 1 << " of " << nz << "\t\r" << std::flush;
            model_fits.read_channel_image(k, image);

            for (int i = 0; i < nx * ny; i++)
            {
                image[i] += gauss_noise(rng);
            }

            model_fits.write_channel_image(k, image);
        }

        return true;
    }

    bool SpectralSimulations::add_vertical_scan_noise_model(fits &model_fits, float sigma_noise, float alpha)
    {

        // Get a random seed from the clock
        // Example from https://cplusplus.com/reference/random/mersenne_twister_engine/seed/
        typedef std::chrono::high_resolution_clock myclock;
        myclock::time_point beginning = myclock::now();

        int nx = model_fits.get_naxes(0);
        int ny = model_fits.get_naxes(1);
        int nz = model_fits.get_naxes(2);

        std::random_device dev;
        std::mt19937 rng(dev());
        std::normal_distribution<float> gauss_noise(0, sigma_noise);

        // Always the same random simulation
        // rng.seed(1);

        // Obtain a seed from the timer
        myclock::duration d = myclock::now() - beginning;
        unsigned s = d.count();
        rng.seed(s);
        cout << "Using random seed: " << s << endl;

        vector<float> image(nx * ny);
        float old_noise = 0;

        for (int k = 0; k < nz; k++)
        {
            cout << "Working on channel: " << k + 1 << " of " << nz << "\t\r" << std::flush;
            model_fits.read_channel_image(k, image);

            for (int i = 0; i < nx; i++)
            {
                for (int j = 0; j < ny; j++)
                {
                    image[i + j * nx] += alpha * old_noise + gauss_noise(rng);
                    old_noise = image[i + j * nx];
                }
            }

            model_fits.write_channel_image(k, image);
        }

        return true;
    }

    bool SpectralSimulations::add_horizontal_scan_noise_model(fits &model_fits, float sigma_noise, float alpha)
    {

        // Get a random seed from the clock
        // Example from https://cplusplus.com/reference/random/mersenne_twister_engine/seed/
        typedef std::chrono::high_resolution_clock myclock;
        myclock::time_point beginning = myclock::now();

        int nx = model_fits.get_naxes(0);
        int ny = model_fits.get_naxes(1);
        int nz = model_fits.get_naxes(2);

        std::random_device dev;
        std::mt19937 rng(dev());
        std::normal_distribution<float> gauss_noise(0, sigma_noise);

        // Always the same random simulation
        // rng.seed(1);

        // Obtain a seed from the timer
        myclock::duration d = myclock::now() - beginning;
        unsigned s = d.count();
        rng.seed(s);
        cout << "Using random seed: " << s << endl;

        vector<float> image(nx * ny);
        float old_noise = 0;

        for (int k = 0; k < nz; k++)
        {
            cout << "Working on channel: " << k + 1 << " of " << nz << "\t\r" << std::flush;
            model_fits.read_channel_image(k, image);

            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    image[i + j * nx] += alpha * old_noise + gauss_noise(rng);
                    old_noise = image[i + j * nx];
                }
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

    bool SpectralSimulations::add_RFI(fits &model_fits, const string &RFI_infile)
    {

        // Open RFI infile
        ifstream in(RFI_infile);
        if (!in.good())
        {
            cout << "Cannot open RFI infile" << endl;
            return false;
        }

        // Read rfi
        vector<rfi> vrfi;
        while (1)
        {
            string s;
            getline(in, s);
            if (s.length() > 0 && s[0] != '#')
            {
                int xp1, xp2, yp1, yp2, ch1, ch2;
                float i;
                stringstream ss(s);
                ss >> xp1;
                ss >> yp1;
                ss >> xp2;
                ss >> yp2;
                ss >> i;
                ss >> ch1;
                ss >> ch2;
                cout << xp1 << " " << xp2 << " " << yp1 << " " << yp2 << " " << i << " " << ch1 << " " << ch2 << endl;
                vrfi.push_back(rfi(xp1, yp1, xp2, yp2, i, ch1, ch2));
            }
            if (in.eof())
                break;
        }

        cout << "Found " << vrfi.size() << " rfi" << endl;

        if (vrfi.size() == 0)
        {
            return true;
        }

        int nx = model_fits.get_naxes(0);
        int ny = model_fits.get_naxes(1);
        int nz = model_fits.get_naxes(2);

        vector<float> image(nx * ny);

        for (int k = 0; k < nz; k++)
        {
            cout << "Working on channel: " << k + 1 << " of " << nz << "\t\r" << std::flush;
            model_fits.read_channel_image(k, image);

            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    int idx = i + nx * j;
                    for (int n = 0; n < vrfi.size(); n++)
                    {
                        image[idx] += vrfi[n].get_intensity(i + 1, j + 1, k + 1);
                    }
                }
            }
            model_fits.write_channel_image(k, image);
        }

        return true;
    }

    bool SpectralSimulations::add_spectral_line(const string &infile, int xp, int yp, int channel, float peak_intensity, float width)
    {

        int status = 0;
        spacew::fits infits;

        if (!infits.open(infile))
        {
            return false;
        }

        if (!add_spectral_line_model(infits, xp, yp, channel, peak_intensity, width))
        {
            return false;
        }

        return true;
    }

}
