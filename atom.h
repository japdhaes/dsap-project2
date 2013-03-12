#ifndef ATOM_H
#define ATOM_H

#define Xunit 3.405

#include <armadillo>
#include <iostream>
#include "cell.h"

using namespace std;
using namespace arma;

class Cell;

class Atom
{
    public:
        Atom(){}
        Atom(vec r, vec v);
        Atom(string type, vec3 r, vec v);
        vec3 getPosition();
        vec3 getVelocity();
        void setPosition(const vec3 &newPosition);
        void setVelocity(const vec3 &newVelocity);
        vec3 getAcceler();
        void setAcceler(const vec3 &newAcceler);

        friend ostream& operator<<( ostream&,  Atom&);

        string chemelement;
        double m;
        int number;

        //used in Class cell
        Atom* previousAtom;
        Atom* nextAtom;
        Cell* currentcell;
        vec3 realposition;
        vec3 initialposition;
    protected:
        vec3 position;
        vec3 velocity;
        vec3 acceler;

};

inline vec3 Atom::getPosition(){
    return position;
}

inline vec3 Atom::getVelocity(){
    return velocity;
}

#endif // ATOM_H
