#include <iostream>
#include "SpectralWeighting.h"

using namespace std;
using namespace spacew;

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cout << "Usage: SpectralWeighting input.txt"<< endl;
        return EXIT_FAILURE;
    }

    SpectralWeighting sp;

    if(!sp.open(argv[1])){
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}