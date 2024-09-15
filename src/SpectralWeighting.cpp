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
        float_nan = std::nanf("nan");
        exponent = 2;
        cout << "Inizialization" << endl;
    }

    SpectralWeighting::~SpectralWeighting()
    {
        cout << "Clean up" << endl;
    }

    bool SpectralWeighting::set_exponent(float e)
    {
        exponent = e;
        return true;
    }

    bool SpectralWeighting::splat(const string &infile, const string &outfile, int bchan, int echan)
    {
        // Open input file and create output file
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
        // Select channel range for splat
        if (echan == 0)
        {
            echan = infits.get_naxes(2);
        }
        else if (echan > 0)
        {
            echan = echan - 1;
        }
        if (bchan > 0)
        {
            bchan = bchan - 1;
        }
        cout << "Splat cube from channel " << bchan + 1 << " to channel " << echan + 1 << endl;

        int nx = infits.get_naxes(0);
        int ny = infits.get_naxes(1);
        int nz = infits.get_naxes(2);
        vector<float> image(nx * ny);
        vector<float> sum_image(nx * ny, 0.0);
        vector<float> wsum_image(nx * ny, 0.0);
        vector<float> splat_image(nx * ny, float_nan);

        for (int k = bchan; k < echan; k++)
        {
            // Loop over channels
            cout << "Working on channel: " << k + 1 << " of " << echan + 1 << "\t\r" << std::flush;

            infits.read_channel_image(k, image);

            for (int idx = 0; idx < nx * ny; idx++)
            {

                if (std::isfinite(image[idx]))
                {
                    sum_image[idx] += image[idx];
                    wsum_image[idx] += 1.0;
                }
            }

            for (int idx = 0; idx < nx * ny; idx++)
            {

                if (std::isfinite(sum_image[idx]) && std::isfinite(wsum_image[idx]) && wsum_image[idx] > 0)
                {
                    splat_image[idx] = sum_image[idx] / wsum_image[idx];
                }
            }
        }
        cout << endl;

        long pix[4];
        pix[0] = 1;
        pix[1] = 1;
        pix[2] = 1;
        pix[3] = 1;
        long npixels = nx * ny;

        // Write the final spat image
        fits_write_pix(outfits.get_fptr(), TFLOAT, pix, npixels, &splat_image[0], &status);

        // Close files
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

    bool SpectralWeighting::gaussian_smoothing(const string &infile, const string &outfile, float sigma)
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

        // Open input file and create output file
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

        cout << "Cloning header" << endl;
        outfits.clone_header(infits);
        outfits.fill(0);

        // Smoothing main loop
        int nx = infits.get_naxes(0);
        int ny = infits.get_naxes(1);
        int nz = infits.get_naxes(2);

        vector<float> image(nx * ny);
        vector<float> smooth_image(nx * ny);

        for (int k = 0; k < nz; k++)
        {
            std::fill(smooth_image.begin(), smooth_image.end(), 0.0);
            cout << "Working on channel: " << k + 1 << " of " << nz << "\t\r" << std::flush;
            infits.read_channel_image(k, image);

            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    for (int ii = 0; ii < m; ii++)
                    {
                        for (int jj = 0; jj < n; jj++)
                        {
                            int xpix = i + ii - int(m / 2.0 + 0.5);
                            int ypix = j + jj - int(n / 2.0 + 0.5);

                            float value = image[xpix + nx * ypix];

                            // If inside image
                            if (std::isfinite(value) && xpix >= 0 && xpix < nx && ypix >= 0 && ypix < ny)
                            {
                                smooth_image[i + nx * j] = smooth_image[i + nx * j] + value * kernel[ii][jj];
                            }
                        }
                    }
                }
            }

            outfits.write_channel_image(k, smooth_image);
        }

        cout << endl;

        // Close files
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

    bool SpectralWeighting::local_noise(const string &infile, const string &outfile, int size)
    {

        // Open input file and create output file
        int status = 0;
        spacew::fits infits;

        if (!infits.open(infile))
        {
            return false;
        }

        cout << "Creating file " << outfile << endl;

        spacew::fits outfits;
        if (!outfits.create(outfile))
        {
            return false;
        }

        outfits.clone_header(infits);
        outfits.fill(0);

        // Smoothing main loop
        int nx = infits.get_naxes(0);
        int ny = infits.get_naxes(1);
        int nz = infits.get_naxes(2);

        vector<float> image(nx * ny);
        vector<float> sigma_image(nx * ny);

        int m = 1 + 2 * size; // It's always odd

        for (int k = 0; k < nz; k++)
        {
            cout << "Working on channel: " << k + 1 << " of " << nz << "\t\r" << std::flush;
            infits.read_channel_image(k, image);
            get_plane_sigma_image(image, nx, ny, sigma_image, m);
            outfits.write_channel_image(k, sigma_image);
        }

        cout << endl;

        // Close files
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

    bool SpectralWeighting::get_plane_sigma_image(vector<float> &image, int nx, int ny, vector<float> &sigma_image, int m)
    {

        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {

                vector<float> values;

                // Extract nearby pixel values
                for (int ii = 0; ii < m; ii++)
                {
                    for (int jj = 0; jj < m; jj++)
                    {
                        int xpix = i + ii - int(m / 2.0 + 0.5);
                        int ypix = j + jj - int(m / 2.0 + 0.5);

                        // If inside image
                        if (xpix >= 0 && xpix < nx && ypix >= 0 && ypix < ny)
                        {
                            float value = image[xpix + nx * ypix];
                            if (std::isfinite(value))
                            {
                                values.push_back(value);
                            }
                        }
                    }
                }

                // Calculate statistic of nearby pixel values
                float sum = 0;
                for (auto v : values)
                {
                    sum += v;
                }
                float avg = float_nan;
                sigma_image[i + nx * j] = float_nan;
                if (values.size() > 0)
                {
                    avg = sum / values.size(); // Average
                    sum = 0;
                    for (auto v : values)
                    {
                        sum += pow((avg - v), 2);
                    }
                    sigma_image[i + nx * j] = sqrt(sum / values.size()); // Standard Deviation
                }
            }
        }

        return true;
    }

    bool SpectralWeighting::local_weights(const string &infile, const string &outfile, int size, int bchan, int echan, float sigma)
    {

        // Smooth input image if sigma>0
        string smooth_infile = infile;
        if (sigma > 0)
        {
            if (!gaussian_smoothing(infile, "!smoothed_" + infile, sigma))
            {
                return false;
            }
            smooth_infile = "smoothed_" + infile;
        }

        // Open input file and create output file
        int status = 0;
        spacew::fits infits;

        if (!infits.open(smooth_infile))
        {
            return false;
        }

        cout << "Creating file " << outfile << endl;

        spacew::fits outfits;
        if (!outfits.create(outfile))
        {
            return false;
        }

        outfits.clone_header(infits);
        outfits.fill(0);

        // Smoothing main loop
        int nx = infits.get_naxes(0);
        int ny = infits.get_naxes(1);
        int nz = infits.get_naxes(2);

        vector<float> image(nx * ny);
        vector<float> sigma_image(nx * ny);

        int m = 1 + 2 * size; // It's always odd

        // Select channel range for Weights' cube
        if (echan == 0)
        {
            echan = infits.get_naxes(2);
        }
        else if (echan > 0)
        {
            echan = echan - 1;
        }
        if (bchan > 0)
        {
            bchan = bchan - 1;
        }
        cout << "Weights' cube from channel " << bchan + 1 << " to channel " << echan + 1 << endl;

        vector<float> mini_splat_image(nx * ny, float_nan);
        vector<float> sum_image(nx * ny, 0.0);
        vector<float> wsum_image(nx * ny, 0.0);

        for (int k = bchan; k < echan; k++)
        {
            cout << "Working on channel: " << k + 1 << " of " << echan + 1 << "\t\r" << std::flush;
            std::fill(mini_splat_image.begin(), mini_splat_image.end(), float_nan);
            std::fill(sum_image.begin(), sum_image.end(), 0.0);
            std::fill(wsum_image.begin(), wsum_image.end(), 0.0);

            // We compute a mini splat of the nearby channels to highlight the presence of faint broad band RFI
            int k1 = std::max(bchan, k - m);
            int k2 = std::min(echan, k + m);

            for (int kk = k1; kk < k2; kk++)
            {
                infits.read_channel_image(kk, image);

                for (int idx = 0; idx < nx * ny; idx++)
                {
                    if (std::isfinite(image[idx]))
                    {
                        sum_image[idx] += image[idx];
                        wsum_image[idx] += 1.0;
                    }
                }

                for (int idx = 0; idx < nx * ny; idx++)
                {
                    if (std::isfinite(sum_image[idx]) && std::isfinite(wsum_image[idx]) && wsum_image[idx] > 0)
                    {
                        mini_splat_image[idx] = sum_image[idx] / wsum_image[idx];
                    }
                }
            }

            get_plane_sigma_image(mini_splat_image, nx, ny, sigma_image, m);
            for (int i = 0; i < sigma_image.size(); i++)
            {
                sigma_image[i] = 1 / pow(sigma_image[i], exponent);
            }
            outfits.write_channel_image(k, sigma_image);
        }

        cout << endl;

        // Close files
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

    bool SpectralWeighting::weighted_splat(const string &infile, const string &outfile, int size, int bchan, int echan, float sigma)
    {
        // Open input file and create output file
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

        // Calculate weights cube and open it
        cout << "Calculating local weights" << endl;
        local_weights(infile, "!weights_" + infile, size, bchan, echan, sigma); // Create a fits cube creating weights

        spacew::fits winfits;

        if (!winfits.open("weights_" + infile))
        {
            return false;
        }

        // Weigthed splat main loop
        // Select channel range for splat
        if (echan == 0)
        {
            echan = infits.get_naxes(2);
        }
        else if (echan > 0)
        {
            echan = echan - 1;
        }
        if (bchan > 0)
        {
            bchan = bchan - 1;
        }
        cout << "Splat cube from channel " << bchan + 1 << " to channel " << echan + 1 << endl;

        int nx = infits.get_naxes(0);
        int ny = infits.get_naxes(1);
        int nz = infits.get_naxes(2);
        vector<float> image(nx * ny);
        vector<float> wimage(nx * ny);
        vector<float> sum_image(nx * ny, 0.0);
        vector<float> wsum_image(nx * ny, 0.0);
        vector<float> splat_image(nx * ny, float_nan);

        for (int k = bchan; k < echan; k++)
        {
            // Loop over channels
            cout << "Working on channel: " << k + 1 << " of " << echan + 1 << "\t\r" << std::flush;

            infits.read_channel_image(k, image);
            winfits.read_channel_image(k, wimage);

            for (int idx = 0; idx < nx * ny; idx++)
            {

                if (std::isfinite(wimage[idx]) && std::isfinite(image[idx]))
                {
                    sum_image[idx] += image[idx] * wimage[idx];
                    wsum_image[idx] += wimage[idx];
                }
            }

            for (int idx = 0; idx < nx * ny; idx++)
            {

                if (std::isfinite(sum_image[idx]) && std::isfinite(wsum_image[idx]) && wsum_image[idx] > 0)
                {
                    splat_image[idx] = sum_image[idx] / wsum_image[idx];
                }
            }
        }
        cout << endl;

        long pix[4];
        pix[0] = 1;
        pix[1] = 1;
        pix[2] = 1;
        pix[3] = 1;
        long npixels = nx * ny;

        // Write the final spat image
        fits_write_pix(outfits.get_fptr(), TFLOAT, pix, npixels, &splat_image[0], &status);

        // Close files
        if (!infits.close()) // Check that worked
        {
            return false;
        }

        if (!winfits.close()) // Check that worked
        {
            return false;
        }

        if (!outfits.close()) // Check that worked
        {
            return false;
        }

        return true;
    }

    bool SpectralWeighting::weighted_merge(const string &filelist, const string &outfile, int size, int bchan, int echan, float sigma)
    {

        // If size=0, no weight is applied

        // Open filelist
        ifstream in(filelist);
        if (!in.good())
        {
            cout << "Cannot open filelist" << endl;
            return false;
        }

        // Read files
        vector<string> files;
        while (1)
        {
            string s;
            getline(in, s);
            if (s.length() > 0 && s[0] != '#')
            {
                files.push_back(s);
            }
            if (in.eof())
                break;
        }

        if (files.size() == 0)
        {
            cout << "No input files found in filelist" << endl;
            return false;
        }

        cout << "Processing " << files.size() << " files" << endl;

        // Produce weights cubes
        if (size > 0.0)
        {
            for (int i = 0; i < files.size(); i++)
            {
                if (!local_weights(files[i], "!weights_" + files[i], size, bchan, echan, sigma))
                {
                    return false;
                }
            }
        }

        // Open input and output files
        int status = 0;
        vector<spacew::fits> vinfits(files.size());
        vector<spacew::fits> vwinfits(files.size());

        for (int i = 0; i < files.size(); i++)
        {
            if (!vinfits[i].open(files[i]))
            {
                return false;
            }
            if (size > 0.0 && !vwinfits[i].open("weights_" + files[i]))
            {
                return false;
            }
        }

        spacew::fits outfits;
        if (!outfits.create(outfile))
        {
            return false;
        }

        outfits.clone_header(vinfits[0]);
        outfits.fill(0);

        // Weigthed Merge main loop
        // Select channel range for merge
        if (echan == 0)
        {
            echan = vinfits[0].get_naxes(2);
        }
        else if (echan > 0)
        {
            echan = echan - 1;
        }
        if (bchan > 0)
        {
            bchan = bchan - 1;
        }
        cout << "Merge cube from channel " << bchan + 1 << " to channel " << echan + 1 << endl;

        int nx = vinfits[0].get_naxes(0);
        int ny = vinfits[0].get_naxes(1);
        int nz = vinfits[0].get_naxes(2);

        for (int k = bchan; k < echan; k++)
        {
            // Loop over channels
            cout << "Working on channel: " << k + 1 << " of " << echan + 1 << "\t\r" << std::flush;

            vector<float> sum_image(nx * ny, 0.0);
            vector<float> wsum_image(nx * ny, 0.0);
            vector<float> merge_image(nx * ny, float_nan);

            // Loop over files
            for (int i = 0; i < files.size(); i++)
            {

                vector<float> image(nx * ny, 0.0);

                float winit = 0.0;
                if (size == 0)
                {
                    winit = 1.0;
                }
                vector<float> wimage(nx * ny, winit);

                vinfits[i].read_channel_image(k, image);
                if (size > 0)
                {
                    vwinfits[i].read_channel_image(k, wimage);
                }

                for (int idx = 0; idx < nx * ny; idx++)
                {
                    if (std::isfinite(wimage[idx]) && std::isfinite(image[idx]))
                    {
                        sum_image[idx] += image[idx] * wimage[idx];
                        wsum_image[idx] += wimage[idx];
                    }
                }
            }
            for (int idx = 0; idx < nx * ny; idx++)
            {

                if (std::isfinite(sum_image[idx]) && std::isfinite(wsum_image[idx]) && wsum_image[idx] > 0)
                {
                    merge_image[idx] = sum_image[idx] / wsum_image[idx];
                }
            }

            long pix[4];
            pix[0] = 1;
            pix[1] = 1;
            pix[2] = k + 1;
            pix[3] = 1;
            long npixels = nx * ny;

            // Write the final merge image
            fits_write_pix(outfits.get_fptr(), TFLOAT, pix, npixels, &merge_image[0], &status);
        }
        cout << endl;

        // Close files
        for (int i = 0; i < files.size(); i++)
        {

            if (!vinfits[i].close()) // Check that worked
            {
                return false;
            }

            if (size > 0 && !vwinfits[i].close()) // Check that worked
            {
                return false;
            }
        }

        if (!outfits.close()) // Check that worked
        {
            return false;
        }

        return true;
    }

    bool SpectralWeighting::sum_fits(const string &infile1, const string &infile2, const string &outfile, float f1, float f2)
    {

        int status = 0;
        spacew::fits infits1;
        spacew::fits infits2;

        if (!infits1.open(infile1))
        {
            return false;
        }

        if (!infits2.open(infile2))
        {
            return false;
        }

        spacew::fits outfits;
        if (!outfits.create(outfile))
        {
            return false;
        }

        outfits.clone_header(infits1);
        outfits.fill(0);

        // Add main loop
        int nx = infits1.get_naxes(0);
        int ny = infits1.get_naxes(1);
        int nz = infits1.get_naxes(2);

        vector<float> image1(nx * ny);
        vector<float> image2(nx * ny);
        vector<float> add_image(nx * ny);

        for (int k = 0; k < nz; k++)
        {
            cout << "Working on channel: " << k + 1 << " of " << nz << "\t\r" << std::flush;
            infits1.read_channel_image(k, image1);
            infits2.read_channel_image(k, image2);

            for (int i = 0; i < nx * ny; i++)
            {
                add_image[i] = f1 * image1[i] + f2 * image2[i];
            }
            outfits.write_channel_image(k, add_image);
        }

        cout << endl;

        // Close files
        if (!infits1.close()) // Check that worked
        {
            return false;
        }

        if (!infits2.close()) // Check that worked
        {
            return false;
        }

        if (!outfits.close()) // Check that worked
        {
            return false;
        }

        return true;
    }

    bool SpectralWeighting::flag_channels(const string &infile, const string &outfile, float sigma_threshold, int bchan, int echan)
    {
        // Open input file and create output file
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
        outfits.fill(float_nan);

        vector<ChannelStatistic> vstat;
        get_spectrum(infile, vstat, bchan, echan);

        int nx = infits.get_naxes(0);
        int ny = infits.get_naxes(1);
        int nz = infits.get_naxes(2);
        vector<float> image(nx * ny);

        if (echan == 0)
        {
            echan = infits.get_naxes(2);
        }
        else if (echan > 0)
        {
            echan = echan - 1;
        }
        if (bchan > 0)
        {
            bchan = bchan - 1;
        }

        for (int k = bchan; k < echan; k++)
        {
            // Loop over channels
            cout << "Working on channel: " << k + 1 << " of " << echan + 1 << "\t\r" << std::flush;
            bool save = false;

            for (int i = 0; i < vstat.size(); i++)
            {
                if (vstat[i].channel == k && vstat[i].sigma < sigma_threshold)
                {
                    save = true;
                }
            }

            if (save)
            {
                infits.read_channel_image(k, image);
                long pix[4];
                pix[0] = 1;
                pix[1] = 1;
                pix[2] = k + 1;
                pix[3] = 1;
                long npixels = nx * ny;
                fits_write_pix(outfits.get_fptr(), TFLOAT, pix, npixels, &image[0], &status);
            }
        }

        // Close files
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

    bool SpectralWeighting::get_spectrum(const string &infile, vector<ChannelStatistic> &vstat, int bchan, int echan)
    {

        // Open input file and create output file
        int status = 0;
        spacew::fits infits;

        if (!infits.open(infile))
        {
            return false;
        }

        int nx = infits.get_naxes(0);
        int ny = infits.get_naxes(1);
        int nz = infits.get_naxes(2);

        if (echan == 0)
        {
            echan = infits.get_naxes(2);
        }
        else if (echan > 0)
        {
            echan = echan - 1;
        }
        if (bchan > 0)
        {
            bchan = bchan - 1;
        }

        for (int k = bchan; k < echan; k++)
        {
            // Loop over channels
            cout << "Working on channel: " << k + 1 << " of " << echan + 1 << "\t\r" << std::flush;

            float average = 0;
            float sigma = 0;
            int npixel = 0;
            vector<float> image;

            infits.read_channel_image(k, image);

            for (int idx = 0; idx < nx * ny; idx++)
            {
                float value = image[idx];

                if (std::isfinite(value))
                {
                    average = average + value;
                    npixel++;
                }
            }

            if (npixel > 0)
            {

                average = average / npixel;

                for (int idx = 0; idx < nx * ny; idx++)
                {
                    float value = image[idx];

                    if (std::isfinite(value))
                    {
                        sigma = sigma + pow((value - average), 2);
                    }
                }
                sigma = sqrt(sigma / npixel);
            }
            else
            {
                average = float_nan;
            }

            ChannelStatistic stat;
            stat.channel = k;
            stat.average = average;
            stat.sigma = sigma;
            stat.npix = npixel;

            vstat.push_back(stat);
        }
        cout << endl;

        // Close files
        if (!infits.close()) // Check that worked
        {
            return false;
        }

        return true;
    }
}
