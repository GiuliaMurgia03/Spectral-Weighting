#pragma once
#include "fitsio.h"
#include <string>
#include <cmath>

using namespace std;

namespace spacew
{

    class SpectralWeighting
    {

        float float_nan;
    public:
        SpectralWeighting();
        ~SpectralWeighting();
        bool splat(const string& infile, const string& outfile, int bchan=0, int echan=0);
        bool weighted_splat(const string &infile, const string &outfile, int size, int bchan=0, int echan=0, float sigma=0);
        bool gaussian_smoothing(const string &infile, const string &outfile, float sigma);
        bool local_noise(const string &infile, const string &outfile, int size);
        bool get_plane_sigma_image(vector <float> &image, int nx, int ny, vector <float> &sigma_image, int size);
        bool local_weights(const string &infile, const string &outfile, int size, int bchan=0, int echan=0, float sigma=0);
        bool weighted_merge(const string &filelist, const string &outfile, int size, int bchan=0, int echan=0, float sigma=0);
        
    };

}