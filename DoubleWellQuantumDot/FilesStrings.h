#pragma once
#include <string>
#include <iomanip>
#include <sstream>
#include <iostream>


class FilesStrings
{
public:
    static double extractDouble(std::string text)
    {
        using namespace std;
        const string digits = "0123456789";
        double x = 0.0;

        unsigned ipos = text.find_first_of(digits);
        if (ipos != string::npos) stringstream(text.substr(ipos)) >> x;
        else                        cout << "Improper input!\n";

        return x;
    }

    static int print_matrix(FILE* f, const gsl_matrix* m)
    {
        int status, n = 0;

        for (size_t i = 0; i < m->size1; i++) {
            for (size_t j = 0; j < m->size2; j++) {
                if ((status = fprintf(f, "%g ", gsl_matrix_get(m, i, j))) < 0)
                    return -1;
                n += status;
            }

            if ((status = fprintf(f, "\n")) < 0)
                return -1;
            n += status;
        }

        return n;
    }
};

