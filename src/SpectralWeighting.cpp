#include <iostream>
#include <fstream>
#include "SpectralWeighting.h"

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

    bool SpectralWeighting::test(const string &filename)
    {

        // Open fits file

        fitsfile *fptr;
        int status = 0;
        fits_open_file(&fptr, filename.c_str(), READONLY, &status);

        if (status) // Check that worked
        {
            fits_report_error(stderr, status);
            return false;
        }

        // Get Header information

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

        // Read pixel value

        long fpixel[naxis];
        fpixel[3] = 1;
        long nelements = 1;
        float array[nelements];

        while (1)
        {
            cout << "Enter x, y, c: ";
            cin >> fpixel[0] >> fpixel[1] >> fpixel[2];
            if (fpixel[0] < 0)
                break;
            fits_read_pix(fptr, TFLOAT, fpixel, nelements, NULL, array, NULL, &status);
            cout << "Pixel value: " << array[0] << endl;
        }

        // Create a fits file and write a test image on it

        fitsfile *ofptr;
        string ofilename = "!test.fits";

        status = 0;
        fits_create_file(&ofptr, ofilename.c_str(), &status);

        if (status) // Check that worked
        {
            fits_report_error(stderr, status);
            return false;
        }

        // Copy the header keywords from input file to output file
        fits_copy_header(fptr, ofptr, &status);

        // Switch the number of channels to 1
        long naxis3 = 1;
        fits_update_key(ofptr, TLONG, "NAXIS3", &naxis3, "Splat plane", &status);

        // Inizialize the values in the image with a linear ramp function
        long pix[4];
        pix[2] = 1;
        pix[3] = 1;
        long nvalues = 1;
        float pixvalues[nvalues];

        for (int j = 0; j < naxes[1]; j++)
        {
            cout << "Working on row: " << j + 1 << " of " << naxes[1] << endl;
            for (int i = 0; i < naxes[0]; i++)
            {
                pix[0] = i+1;
                pix[1] = j+1;
                pixvalues[0] = i + j;
                fits_write_pix(ofptr, TFLOAT, pix, nvalues, pixvalues, &status);
            }
        }

        // Close files
        fits_close_file(fptr, &status);
        if (status) // Check that worked
        {
            fits_report_error(stderr, status);
            return false;
        }

        fits_close_file(ofptr, &status);
        if (status) // Check that worked
        {
            fits_report_error(stderr, status);
            return false;
        }

        return true;
    }

}
