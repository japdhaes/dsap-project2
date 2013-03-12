#include "cell.h"


Cell::Cell()
{
    this->first=NULL;
}

void Cell::insertElement(Atom *atom)
{
    atom->currentcell=this;
    //no atoms in this cell so far
    if(this->first==NULL){
        this->first=atom;
        atom->previousAtom=NULL;
        atom->nextAtom=NULL;
        return;
    }

    //already at least 1 atom in the cell
    atom->nextAtom=this->first;
    this->first->previousAtom=atom;
    this->first=atom;
    atom->previousAtom=NULL;


    if(this->isAtomInCell(atom)){
        //debugging << "joepie atom in cell!"<< endl;
    }
    //debugging.close();
}

void Cell::removeelement(Atom *atom){
    Atom *otheratom = first;
    //if first element in the cell is already this atom
    if(otheratom==atom){
        this->first=atom->nextAtom;
        atom->nextAtom=NULL;
        atom->previousAtom=NULL;
        return;
    }

    //else it's somewhere in the middle or at the end
    (atom->previousAtom)->nextAtom=atom->nextAtom;
    if(atom->nextAtom!=NULL){
        (atom->nextAtom)->previousAtom=atom->previousAtom;
    }
    atom->nextAtom=NULL;
    atom->previousAtom=NULL;
}

int Cell::countAtomsInCell(){
    Atom *atom = first;
    int answer=0;
    while(atom!=NULL){
        answer++;
        atom=atom->nextAtom;
    }
    return answer;
}

bool Cell::isAtomInCell(Atom *atom){
    vec3 r = atom->getPosition();
    int intpos[3];

    for(int j=0; j<3; j++){
        intpos[j]=int(r(j)/this->vectorBC(j));
    }
    if(x==intpos[0] && y==intpos[1] && z==intpos[2]){
        return true;
    }
    return false;
}
