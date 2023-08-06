#include <iostream>
#include "fits.h"

using namespace std;

namespace spacew
{

    fits::fits()
    {

        naxis = 0;
        naxes = nullptr;
        fptr = nullptr;
        status = 0;
    }

    fits::~fits()
    {

        if (naxes)
        {
            delete[] naxes;
            naxes = nullptr;
        }
    }

    int fits::get_naxis() const
    {
        return naxis;
    }

    long fits::get_naxes(const int &i) const
    {
        return naxes[i];
    }

    bool fits::set_naxes(const int &i, const long &value)
    {
        if (i >= 0 && i < naxis)
        {
            naxes[i] = value;
            string s = "NAXIS" + std::to_string(i + 1);
            status = 0;
            fits_update_key(fptr, TLONG, s.c_str(), &naxes[i], "Axis size", &status);

            if (status) // Check that worked
            {
                fits_report_error(stderr, status);
                return false;
            }
            return true;
        }

        return false;
    }

    bool fits::open(const string &infile)
    {
        status = 0;
        fits_open_file(&fptr, infile.c_str(), READONLY, &status);

        if (status) // Check that worked
        {
            fits_report_error(stderr, status);
            return false;
        }

        get_header_info();
        return true;
    }

    bool fits::get_header_info()
    {
        // Get Header information
        int hdunum;
        fits_get_num_hdus(fptr, &hdunum, &status);
        cout << "Number of HDU in fits file: " << hdunum << endl;

        int bitpix;
        fits_get_img_type(fptr, &bitpix, &status);
        cout << "BitPix is: " << bitpix << endl;

        fits_get_img_dim(fptr, &naxis, &status);
        cout << "The number of axis is: " << naxis << endl;

        naxes = new long[naxis]; // Array of length naxis
        fits_get_img_size(fptr, naxis, naxes, &status);
        for (int i = 0; i < naxis; i++)
            cout << "Axis " << i << " size is: " << naxes[i] << endl;

        return true;
    }

    bool fits::create(const string &outfile)
    {

        status = 0;
        fits_create_file(&fptr, outfile.c_str(), &status);

        if (status) // Check that worked
        {
            fits_report_error(stderr, status);
            return false;
        }

        return true;
    }

    bool fits::close()
    {

        // Close file
        status = 0;
        fits_close_file(fptr, &status);
        if (status) // Check that worked
        {
            fits_report_error(stderr, status);
            return false;
        }

        return true;
    }

    bool fits::clone_header(const fits &infits)
    {
        // Copy the header keywords from input file to output file
        status = 0;
        fits_copy_header(infits.get_fptr(), fptr, &status);

        if (status) // Check that worked
        {
            fits_report_error(stderr, status);
            return false;
        }

        naxis = infits.get_naxis();
        if (naxes)
        {
            delete[] naxes;
        }

        naxes = new long[naxis];
        for (int i = 0; i < naxis; i++)
        {
            naxes[i] = infits.get_naxes(i);
        }

        return true;
    }

    fitsfile *fits::get_fptr() const
    {

        return fptr;
    }

    bool fits::fill(const float &value)
    {
        // Inizialize the values in the image with a value provided by the user
        long pix[naxis];
        for (int i = 0; i < naxis; i++)
        {
            pix[i] = 1;
        }

        for (int i = 0; i < naxis; i++)
        {
            cout << "Naxes " << i << " : " << naxes[i] << endl;
        }

        long nvalues = naxes[0] * naxes[1];

        float pixvalues[nvalues];

        for (int i = 0; i < nvalues; i++)
        {
            pixvalues[i] = value;
        }

        if (naxis >= 3)
        {
            for (int c = 0; c < naxes[2]; c++)
            {
                pix[2] = c+1;
                fits_write_pix(fptr, TFLOAT, pix, nvalues, pixvalues, &status);
            }
        }
        else if (naxis == 2)
        {
            fits_write_pix(fptr, TFLOAT, pix, nvalues, pixvalues, &status);
        }

        if (status) // Check that worked
        {
            fits_report_error(stderr, status);
            return false;
        }
        return true;
    }

    bool fits::read_channel_image(const int &channel, vector<float> &image)
    {

        long pix[naxis];
        status = 0;

        for (int i = 0; i < naxis; i++)
        {
            pix[i] = 1;
        }

        pix[0] = 1;
        pix[1] = 1;
        pix[2] = channel + 1;
        long nvalues = naxes[0] * naxes[1];

        fits_read_pix(fptr, TFLOAT, pix, nvalues, NULL, &image[0], NULL, &status);

        if (status) // Check that worked
        {
            fits_report_error(stderr, status);
            return false;
        }

        return true;
    }

    bool fits::write_channel_image(const int &channel, vector<float> &image)
    {

        long pix[naxis];
        status = 0;

        for (int i = 0; i < naxis; i++)
        {
            pix[i] = 1;
        }

        pix[0] = 1;
        pix[1] = 1;
        pix[2] = channel + 1;
        long nvalues = naxes[0] * naxes[1];

        fits_write_pix(fptr, TFLOAT, pix, nvalues, &image[0], &status);

        if (status) // Check that worked
        {
            fits_report_error(stderr, status);
            return false;
        }

        return true;
    }

}