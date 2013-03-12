#ifndef VERLETALGO2_H
#define VERLETALGO2_H

#include <armadillo>
#include "crystal.h"

const double cutoffacceleration=200;
//time unit of h =2.1569 *  10^3 fs

using namespace arma;
class VerletAlgo2
{
public:
    ofstream debugging;
    VerletAlgo2(Crystal *crystal, double _h);
    void integrate(bool);
    void integrate_noapprox();
    void integrateAtom(Atom *atom, vec3 boundvec);

    double h;
    Crystal *crystall;

    void updatePosition(Atom *atom);
    void updateVelocity(Atom *atom);
    void updateAcceler(Atom *atom);
    void updateAccelerNoApprox(Atom *atom);

    void calcForce(Atom *atom, Atom *otheratom);
    void findXYZCellIndices(int *nrXYZ, int *nrX, int *nrY, int *nrZ);
    vec3 findClosestPosition(vec3 &position, vec3 &otherposition);
    vec3 boundCheck(vec3 &position);
    double LJpotential(vec3 relvec);
    void thermostatAnders();
    void thermostatBerendsen();
    void thermostatAndersen();
};

#endif // VERLETALGO2_H

