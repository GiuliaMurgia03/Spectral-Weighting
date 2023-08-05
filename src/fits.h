#pragma once
#include "fitsio.h"
#include <string>

using namespace std;

namespace spacew
{

    class fits
    {

    private:
        fitsfile *fptr;
        int naxis;
        long *naxes;
        int status;

    public:
        fits();
        ~fits();
        bool open(const string &infile);
        bool get_header_info();
        fitsfile* get_fptr() const;
        bool create(const string &outfile);
        bool close();
        bool clone_header(const fits &infits);
        bool fill(const float &value);
        int get_naxis() const;     // Returns the number of axis
        long get_naxes(const int& i) const; // Returns the size of each axis
        bool set_naxes(const int& i,const long& value);
    };

}