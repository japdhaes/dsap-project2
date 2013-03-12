#ifndef CRYSTAL_H
#define CRYSTAL_H

#include <iostream>
#include <vector>
#include "atom.h"
#include "cell.h"
#include "zignor.h"
#include "zigrandom.h"
#include "lib.h"

using namespace std;

//3.405 Angstrom
const double xunit  =3.405;
const double xunitSI=3.405e-10;
const double tunitSI=2.1569e-12;
const double vunitSI=xunitSI/tunitSI;
//119.74 kelvin
const double tempunit =119.74;
//0.010318 eV
const double energyunit=0.010318;

//in SI units = Pa
const double pressureunit=energyunit*1.6e-19/(xunit*xunit*xunit*1e-30);


class Crystal
{
    public:
        //ofstream* debugging2;

        Crystal()
        {
        }
        Crystal(unsigned int nc, double _b, int& seed, double _temperature);
        Crystal(string &filename);

        vector<Atom*> allatoms;
        vector<vector<vector<Cell> > > allcells;

        int numberofatoms;
        int nc;
        double beginenergy;
        int counter;
        vec3 boundary;

        double pressure;
        double volume;
        double density;

        int countAtoms();
        friend ostream& operator<<( ostream&, const Crystal&);

        vec3 vectorBC;
        double b;
        double energy;
        double pe;
        double ke;
        double inittemp;
        //mean square displacement
        double msqdplm;
        void setvectorBC(double desiredwidth);
        void initializeAtoms(double _temperature);
        void addAllAtomsToCells();
        void initializeCells();
        void findCellOfAtom(Atom *atom, int &x, int &y, int &z);

        double temperature();
        void removeCrystalMomentum();
        vec3 findClosestPosition(vec3 position, vec3 otherposition);
        void radialDistFunction();

        void createCylinderPore(double r);
        void boundCheck(Atom *atom);
        long idum;
        void createSphere(double begin, double end);
};

#endif // CRYSTAL_H
