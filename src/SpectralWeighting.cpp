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

    bool SpectralWeighting::open(const string &filename)
    {
        ifstream in;
        in.open(filename);

        if (!in.good())
        {
            cout << "Cannot open " << filename << endl;
        }

        fitsfile *fptr;
        int status = 0;
        fits_open_file(&fptr, filename.c_str(), READONLY, &status);

        if (status)
        {
            fits_report_error(stderr, status);
            return false;
        }

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

        long fpixel[naxis];
        fpixel[3] = 1;
        long nelements = 1;
        float array[nelements];

        while (1)
        {
            cout << "Enter x, y, c: ";
            cin >> fpixel[0] >> fpixel[1] >> fpixel[2];
            fits_read_pix(fptr, TFLOAT, fpixel, nelements, NULL, array, NULL, &status);
            cout << "Pixel value: " << array[0] << endl;
        }

        return true;
    }

}
