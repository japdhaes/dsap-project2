#ifndef PRINTING_H
#define PRINTING_H

#include "crystal.h"
#include <sstream>
#include <iostream>
#include <iomanip>

using namespace std;

class Printing
{
public:
    Printing();
    void printing(Crystal &crystal);
    void printvelocities(Crystal &crystal, int j);
    string createname(int number);
};

#endif // PRINTING_H
