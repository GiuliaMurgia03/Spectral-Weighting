#include <iostream>
#include <fstream>
#include <cmath>
#include "SpectralWeighting.h"
#include "fits.h"

using namespace std;

namespace spacew
{

    SpectralWeighting::SpectralWeighting()
    {
        cout << "Inizialization" << endl;
    }

    SpectralWeighting::~SpectralWeighting()
    {
        cout << "Clean up" << endl;
    }

    bool SpectralWeighting::splat(const string &infile, const string &outfile)
    {
        // Open input file and create output file
        int status = 0;
        spacew::fits infits;

        if (!infits.open(infile))
        {
            return false;
        }

        spacew::fits outfits;
        if (!outfits.create(outfile))
        {
            return false;
        }

        outfits.clone_header(infits);
        outfits.set_naxes(2, 1);
        outfits.fill(0);

        // Splat main loop
        long pix[4];
        pix[0] = 1;
        pix[1] = 1;
        pix[2] = 1;
        pix[3] = 1;
        float average[1];
        float channelvalue[1];

        for (int j = 0; j < infits.get_naxes(1); j++)
        {
            cout << "Working on row: " << j + 1 << " of " << infits.get_naxes(1) << endl;
            for (int i = 0; i < infits.get_naxes(0); i++)
            {
                pix[0] = i + 1;
                pix[1] = j + 1;

                average[0] = 0;
                float wsum = 0;
                float sum = 0;

                for (int k = 0; k < infits.get_naxes(2); k++)
                { // Loop over channels
                    pix[2] = k + 1;
                    fits_read_pix(infits.get_fptr(), TFLOAT, pix, 1, NULL, channelvalue, NULL, &status);

                    float w = 1;

                    if (std::isfinite(channelvalue[0]))
                    {
                        sum = sum + w * channelvalue[0];
                        wsum = wsum + w;
                    }
                }

                if (wsum > 0)
                {
                    average[0] = sum / wsum;
                }
                else
                {
                    average[0] = NULL;
                }

                pix[2] = 1; // Reset channel coordinate to 1

                fits_write_pix(outfits.get_fptr(), TFLOAT, pix, 1, average, &status);
            }
        }

        // Close files
        if (!infits.close()) // Check that worked
        {
            return false;
        }

        if (!outfits.close()) // Check that worked
        {
            return false;
        }

        return true;
    }

    bool SpectralWeighting::gaussian_smoothing(const string &infile, const string &outfile, float sigma)
    {

        // Compute the smoothing Kernel
        int m = 1 + 2 * int(3.0 * sigma); // It's always odd
        int n = 1 + 2 * int(3.0 * sigma);

        float kernel[m][n];
        float xc = m / 2.0 + 0.5;
        float yc = n / 2.0 + 0.5;

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                float xp = i + 1;
                float yp = j + 1;
                float d2 = pow(xp - xc, 2) + pow(yp - yc, 2);
                kernel[i][j] = exp(-0.5 * (d2 / pow(sigma, 2)));
            }
        }

        // Open input file and create output file
        int status = 0;
        spacew::fits infits;

        if (!infits.open(infile))
        {
            return false;
        }

        cout << "Creating file " << outfile << endl;

        spacew::fits outfits;
        if (!outfits.create(outfile))
        {
            return false;
        }

        outfits.clone_header(infits);
        outfits.fill(0);

        // Smoothing main loop
        int nx = infits.get_naxes(0);
        int ny = infits.get_naxes(1);
        int nz = infits.get_naxes(2);

        vector<float> image(nx * ny);
        vector<float> smooth_image(nx * ny);

        for (int k = 0; k < nz; k++)
        {
            cout << "Working on channel: " << k + 1 << " of " << nz << endl;
            infits.read_channel_image(k, image);

            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    for (int ii = 0; ii < m; ii++)
                    {
                        for (int jj = 0; jj < n; jj++)
                        {
                            int xpix = i + ii - int(m / 2.0 + 0.5);
                            int ypix = j + jj - int(n / 2.0 + 0.5);

                            // If inside image
                            if (xpix >= 0 && xpix < nx && ypix >= 0 && ypix < ny)
                            {
                                smooth_image[i + nx * j] = smooth_image[i + nx * j] + image[xpix + nx * ypix] * kernel[ii][jj];
                            }
                        }
                    }
                }
            }

            outfits.write_channel_image(k, smooth_image);
        }

        // Close files
        if (!infits.close()) // Check that worked
        {
            return false;
        }

        if (!outfits.close()) // Check that worked
        {
            return false;
        }

        return true;
    }

    bool SpectralWeighting::local_noise(const string &infile, const string &outfile, int size)
    {

        // Open input file and create output file
        int status = 0;
        spacew::fits infits;

        if (!infits.open(infile))
        {
            return false;
        }

        cout << "Creating file " << outfile << endl;

        spacew::fits outfits;
        if (!outfits.create(outfile))
        {
            return false;
        }

        outfits.clone_header(infits);
        outfits.fill(0);

        // Smoothing main loop
        int nx = infits.get_naxes(0);
        int ny = infits.get_naxes(1);
        int nz = infits.get_naxes(2);

        vector<float> image(nx * ny);
        vector<float> sigma_image(nx * ny);

        int m = 1 + 2 * size; // It's always odd

        for (int k = 0; k < nz; k++)
        {
            cout << "Working on channel: " << k + 1 << " of " << nz << endl;
            infits.read_channel_image(k, image);
            get_plane_sigma_image(image, nx, ny, sigma_image, m);
            outfits.write_channel_image(k, sigma_image);
        }

        // Close files
        if (!infits.close()) // Check that worked
        {
            return false;
        }

        if (!outfits.close()) // Check that worked
        {
            return false;
        }

        return true;
    }

    bool SpectralWeighting::get_plane_sigma_image(vector<float> &image, int nx, int ny, vector<float> &sigma_image, int m)
    {

        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {

                vector<float> values;

                // Extract nearby pixel values
                for (int ii = 0; ii < m; ii++)
                {
                    for (int jj = 0; jj < m; jj++)
                    {
                        int xpix = i + ii - int(m / 2.0 + 0.5);
                        int ypix = j + jj - int(m / 2.0 + 0.5);

                        // If inside image
                        if (xpix >= 0 && xpix < nx && ypix >= 0 && ypix < ny)
                        {
                            values.push_back(image[xpix + nx * ypix]);
                        }
                    }
                }

                // Calculate statistic of nearby pixel values
                float sum = 0;
                for (auto v : values)
                {
                    sum += v;
                }
                float avg = NULL;
                sigma_image[i + nx * j] = NULL;
                if (values.size() > 0)
                {
                    avg = sum / values.size(); // Average
                    sum = 0;
                    for (auto v : values)
                    {
                        sum += pow((avg - v), 2);
                    }
                    sigma_image[i + nx * j] = sqrt(sum / values.size()); // Standard Deviation
                }
            }
        }

        return true;
    }

    bool SpectralWeighting::local_weights(const string &infile, const string &outfile, int size)
    {

        // Open input file and create output file
        int status = 0;
        spacew::fits infits;

        if (!infits.open(infile))
        {
            return false;
        }

        cout << "Creating file " << outfile << endl;

        spacew::fits outfits;
        if (!outfits.create(outfile))
        {
            return false;
        }

        outfits.clone_header(infits);
        outfits.fill(0);

        // Smoothing main loop
        int nx = infits.get_naxes(0);
        int ny = infits.get_naxes(1);
        int nz = infits.get_naxes(2);

        vector<float> image(nx * ny);
        vector<float> sigma_image(nx * ny);

        int m = 1 + 2 * size; // It's always odd

        for (int k = 0; k < nz; k++)
        {
            cout << "Working on channel: " << k + 1 << " of " << nz << endl;
            infits.read_channel_image(k, image);
            get_plane_sigma_image(image, nx, ny, sigma_image, m);
            for (int i = 0; i < sigma_image.size(); i++)
            {
                sigma_image[i] = 1 / pow(sigma_image[i], 2); 
            }
            outfits.write_channel_image(k, sigma_image);
        }

        // Close files
        if (!infits.close()) // Check that worked
        {
            return false;
        }

        if (!outfits.close()) // Check that worked
        {
            return false;
        }

        return true;
    }



    bool SpectralWeighting::weighted_splat(const string &infile, const string &outfile, int size)
    {
        // Open input file and create output file
        int status = 0;
        spacew::fits infits;

        if (!infits.open(infile))
        {
            return false;
        }

        spacew::fits outfits;
        if (!outfits.create(outfile))
        {
            return false;
        }

        outfits.clone_header(infits);
        outfits.set_naxes(2, 1);
        outfits.fill(0);

        // Calculate weights cube and open it 
        local_weights(infile, "weights_"+infile, size);

        spacew::fits winfits;

        if (!winfits.open("weights_"+infile))
        {
            return false;
        }

        // Weigthed splat main loop
        long pix[4];
        pix[0] = 1;
        pix[1] = 1;
        pix[2] = 1;
        pix[3] = 1;
        float average[1];
        float channelvalue[1];
        float wchannelvalue[1];

        for (int j = 0; j < infits.get_naxes(1); j++)
        {
            cout << "Working on row: " << j + 1 << " of " << infits.get_naxes(1) << endl;
            for (int i = 0; i < infits.get_naxes(0); i++)
            {
                pix[0] = i + 1;
                pix[1] = j + 1;

                average[0] = 0;
                float wsum = 0;
                float sum = 0;

                for (int k = 0; k < infits.get_naxes(2); k++)
                { // Loop over channels
                    pix[2] = k + 1;
                    // Read actual channel value
                    fits_read_pix(infits.get_fptr(), TFLOAT, pix, 1, NULL, channelvalue, NULL, &status);

                    // Read weight of channel value
                    fits_read_pix(winfits.get_fptr(), TFLOAT, pix, 1, NULL, wchannelvalue, NULL, &status);

                    if (std::isfinite(channelvalue[0]) && std::isfinite(wchannelvalue[0]))
                    {
                        sum = sum + wchannelvalue[0] * channelvalue[0];
                        wsum = wsum + wchannelvalue[0];
                    }
                }

                if (wsum > 0)
                {
                    average[0] = sum / wsum;
                }
                else
                {
                    average[0] = NULL;
                }

                pix[2] = 1; // Reset channel coordinate to 1

                fits_write_pix(outfits.get_fptr(), TFLOAT, pix, 1, average, &status);
            }
        }

        // Close files
        if (!infits.close()) // Check that worked
        {
            return false;
        }

        if (!winfits.close()) // Check that worked
        {
            return false;
        }

        if (!outfits.close()) // Check that worked
        {
            return false;
        }

        return true;
    }

}
