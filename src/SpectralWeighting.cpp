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
            fits_report_error(stderr, status);

        return true;
    }

}
