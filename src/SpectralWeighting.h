#pragma once
#include "fitsio.h"
#include <string>

using namespace std;

namespace spacew
{

    class SpectralWeighting
    {
    public:
        SpectralWeighting();
        ~SpectralWeighting();
        bool splat(const string& infile, const string& outfile);
        bool gaussian_smoothing(const string &infile, const string &outfile, float sigma);
        bool local_noise(const string &infile, const string &outfile, int size);
        bool get_plane_sigma_image(vector <float> &image, int nx, int ny, vector <float> &sigma_image, int size);
        bool local_weights(const string &infile, const string &outfile, int size);
        bool weighted_splat(const string &infile, const string &outfile, int size);
    };

}