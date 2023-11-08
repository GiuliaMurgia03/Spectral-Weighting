#pragma once
#include "fitsio.h"
#include <string>
#include <cmath>
#include <random>
#include <vector>
#include "fits.h"
#include "rfi.h"

using namespace std;

namespace spacew
{

    class SpectralSimulations
    {
    private:
        vector<float> flux_bins;
        float min_flux;
        float max_flux;
        float ds;
        string noise_model;

    public:
        SpectralSimulations();
        ~SpectralSimulations();
        bool set_noise_model(const string &model);
        bool init_flux_bins(const float &s_min, const float &s_max, const float &delta_s);
        bool create_outfits(const string &outfile, int nx, int ny, int nz, int nsources, float sigma_noise, string rfi_infile = "");
        bool add_pointsources_model(fits &model_fits, int nsources, double spectral_index = 0);
        bool add_spectral_line_model(fits &model_fits, int xp, int yp, int channel, float peak_intensity, float width);
        bool add_noise_model(fits &model_fits, float sigma_noise);
        bool add_vertical_scan_noise_model(fits &model_fits, float sigma_noise, float alpha = 0.8);
        bool add_horizontal_scan_noise_model(fits &model_fits, float sigma_noise, float alpha = 0.8);
        bool add_RFI(fits &model_fits, const string &RFI_infile);
        float extract_random_flux();
        bool add_spectral_line(const string &infile, int xp, int yp, int channel, float peak_intensity, float width);
    };

}