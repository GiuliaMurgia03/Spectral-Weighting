#pragma once
#include "fitsio.h"
#include <string>
#include <cmath>
#include <random>
#include "fits.h"

using namespace std;

namespace spacew
{

    class SpectralSimulations
    {
        private:
        vector <float> flux_bins;
        float min_flux;
        float max_flux;
        float ds;
        public:
        SpectralSimulations();
        ~SpectralSimulations();
        bool init_flux_bins(const float& s_min, const float& s_max, const float& delta_s);
        bool create_outfits(const string &outfile, int nx, int ny, int nz, int nsources, float sigma_noise);
        bool add_pointsources_model(fits &model_fits, int nsources, float sigma_noise);
        float extract_random_flux();


      
        
    };

}