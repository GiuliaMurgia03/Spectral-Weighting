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

                    float w = 1; // ADD ROUTINE TO CALCULATE WEIGHT

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

    bool SpectralWeighting::gaussian_smoothing(fitsfile *fptr, const string &outfile, float sigma)
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

        // Get Header information
        int status = 0;
        int hdunum;
        fits_get_num_hdus(fptr, &hdunum, &status);
        cout << "Number of HDU in fits file: " << hdunum << endl;

        int bitpix;
        fits_get_img_type(fptr, &bitpix, &status);
        cout << "BitPix is: " << bitpix << endl;

        int naxis;
        fits_get_img_dim(fptr, &naxis, &status);
        cout << "The number of axis is: " << naxis << endl;

        long naxes[naxis]; // Array
        fits_get_img_size(fptr, naxis, naxes, &status);
        for (int i = 0; i < naxis; i++)
            cout << "Axis " << i << " size is: " << naxes[i] << endl;

        // Create a fits file
        fitsfile *ofptr;
        status = 0;
        fits_create_file(&ofptr, outfile.c_str(), &status);

        if (status) // Check that worked
        {
            fits_report_error(stderr, status);
            return false;
        }

        // Copy the header keywords from input file to output file
        fits_copy_header(fptr, ofptr, &status);

        // Smoothing main loop
        long pix[4];
        pix[0] = 1;
        pix[1] = 1;
        pix[2] = 1;
        pix[3] = 1;

        long smooth_pix[4];
        smooth_pix[0] = 1;
        smooth_pix[1] = 1;
        smooth_pix[2] = 1;
        smooth_pix[3] = 1;
        float smooth[1];
        float channelvalue[1];

        // Inizialize all pixels in all channels to 0
        for (int k = 0; k < naxes[2]; k++)
        {
            for (int j = 0; j < naxes[1]; j++)
            {
                for (int i = 0; i < naxes[0]; i++)
                {

                    pix[0] = i + 1;
                    pix[1] = j + 1;
                    pix[2] = k + 1;
                    channelvalue[0] = 0;

                    fits_write_pix(ofptr, TFLOAT, pix, 1, channelvalue, &status);
                }
            }
        }

        for (int k = 0; k < naxes[2]; k++)
        {
            cout << "Working on channel: " << k + 1 << " of " << naxes[2] << endl;
            for (int j = 0; j < naxes[1]; j++)
            {
                for (int i = 0; i < naxes[0]; i++)
                {
                    pix[0] = i + 1;
                    pix[1] = j + 1;
                    pix[2] = k + 1;

                    smooth_pix[2] = k + 1;

                    fits_read_pix(ofptr, TFLOAT, pix, 1, NULL, channelvalue, NULL, &status);

                    smooth[0] = channelvalue[0];

                    for (int ii = 0; ii < m; ii++)
                    {
                        for (int jj = 0; jj < n; jj++)
                        {
                            smooth_pix[0] = i + ii - int(m / 2.0 + 0.5) + 1;
                            smooth_pix[1] = j + jj - int(n / 2.0 + 0.5) + 1;
                            // If inside image
                            if (smooth_pix[0] > 0 && smooth_pix[0] <= naxes[0] && smooth_pix[1] > 0 && smooth_pix[1] <= naxes[1])
                            {

                                fits_read_pix(fptr, TFLOAT, smooth_pix, 1, NULL, channelvalue, NULL, &status);
                                smooth[0] = smooth[0] + channelvalue[0] * kernel[ii][jj];
                            }
                        }
                    }

                    fits_write_pix(ofptr, TFLOAT, pix, 1, smooth, &status);
                }
            }
        }

        // Close output
        fits_close_file(ofptr, &status);
        if (status) // Check that worked
        {
            fits_report_error(stderr, status);
            return false;
        }

        return true;
    }

}
