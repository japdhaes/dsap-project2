#ifndef CELL_H
#define CELL_H

#include "atom.h"
#include <armadillo>
class Atom;

using namespace arma;

class Cell
{
public:
    Cell();
    Cell(int _x, int _y, int _z):x(_x), y(_y), z(_z){}
    Atom* first;
    double visited;

    void insertElement(Atom* atom);
    void removeelement(Atom *atom);
    void findCell(Atom* atom, double &x, double &y, double &z);
    void testing();

    //position of left corner of the cell
    int x, y, z;
    vec3 vectorBC;
    bool isAtomInCell(Atom* atom);
    int countAtomsInCell();
    int countAtoms();
};

#endif // CELL_H
