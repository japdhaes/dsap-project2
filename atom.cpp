#include "atom.h"

Atom::Atom(vec position, vec v): position(position), velocity(v){
    this->chemelement="Ar";
    vec3 r; r.fill(0);
    this->setAcceler(r);

    this->previousAtom=NULL;
    this->nextAtom=NULL;
    this->m=39.948*1.66e-27;
    this->realposition=position;
    this->initialposition=position;
}

Atom::Atom(string type, vec3 r, vec v){
    this->chemelement=type;
    this->position=r;
    this->velocity=v;
    vec3 nullvector; nullvector.zeros();
    this->setAcceler(nullvector);
    this->previousAtom=NULL;
    this->nextAtom=NULL;
    this->m=39.948*1.66e-27;
    this->realposition=r;
    this->initialposition=r;
}



void Atom::setPosition(const vec3 &newPosition){
    this->position=newPosition;
}

void Atom::setVelocity(const vec3 &newVelocity){
    this->velocity=newVelocity;
}

vec3 Atom::getAcceler()
{
    return this->acceler;
}

void Atom::setAcceler(const vec3 &newAcceler)
{
    this->acceler=newAcceler;
}

ostream& operator<< (ostream& os , Atom& atom){
    //first write down the total number of atoms in the simulated crystal
    vec3 position=atom.getPosition();
    vec3 velocity=atom.getVelocity();
    //os << position(0)*xunit<< " " << position(1)*xunit << " " << position(2)*xunit << " " << velocity(0)<< " " << velocity(1)<< " "<< velocity(2)<< endl;
    os << position(0)<< " " << position(1) << " " << position(2) << " " << velocity(0)<< " " << velocity(1)<< " "<< velocity(2)<< endl;

    return os;
}
