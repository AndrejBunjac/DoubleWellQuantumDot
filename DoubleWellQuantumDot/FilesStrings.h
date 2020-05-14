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
};

