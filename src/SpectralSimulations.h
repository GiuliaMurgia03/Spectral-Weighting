#pragma once
#include "fitsio.h"
#include <string>
#include <cmath>
#include <random>
#include "fits.h"
#include "rfi.h"

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
        bool init_flux_bins(const float& s_min, const float& s_max, const float &delta_s);
        bool create_outfits(const string &outfile, int nx, int ny, int nz, int nsources, float sigma_nois, string rfi_infile="");
        bool add_pointsources_model(fits &model_fits, int nsources);
        bool add_noise_model(fits &model_fits, float sigma_noise);
        bool add_RFI(fits &model_fits, const string &RFI_infile);
        float extract_random_flux();


      
        
    };

}