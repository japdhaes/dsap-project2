#include "verletalgo2.h"

VerletAlgo2::VerletAlgo2(Crystal *crystal, double _h)
{
    this->debugging.open("/home/jonathan/projectsFSAP/project2/project2/debuglog.txt");
    this->crystall=crystal;
    this->h=_h;
}

void VerletAlgo2::integrate(bool thermalize){
    //ofstream debugging;
    //debugging.open("/home/jonathan/projectsFSAP/project1/project1/debuglog2.txt", ios::app);
    crystall->energy=0;
    crystall->pressure=0;

    for(unsigned int i=0; i<crystall->allcells.size(); i++){
        for(unsigned int j=0; j<crystall->allcells.at(i).size();j++){
            for(unsigned int k=0; k<crystall->allcells.at(i).at(j).size(); k++){
                crystall->allcells.at(i).at(j).at(k).visited=false;
            }
        }
    }

    for(unsigned int i=0; i<crystall->allatoms.size(); i++){
        updateVelocity(crystall->allatoms[i]);
        updatePosition(crystall->allatoms[i]);


        vec3 nullvector; nullvector.zeros();
        crystall->allatoms[i]->setAcceler(nullvector);

    }
    this->crystall->msqdplm/=this->crystall->numberofatoms;
//    cout << this->crystall->msqdplm <<endl;
// cout << crystall->temperature()<<endl;
    for(unsigned int i=0; i<crystall->allcells.size(); i++){
        for(unsigned int j=0; j<crystall->allcells.at(i).size();j++){
            for(unsigned int k=0; k<crystall->allcells.at(i).at(j).size(); k++){
                crystall->allcells.at(i).at(j).at(k).visited=true;
                Atom *atom = crystall->allcells.at(i).at(j).at(k).first;
                while(atom!=NULL){
                    updateAcceler(atom);
                    atom=atom->nextAtom;
                }

            }
        }
    }
    crystall->pressure/=3*crystall->volume;
    crystall->pressure+=crystall->density*crystall->temperature();
// crystall->pressure*=pressureunit;
    crystall->pe=crystall->energy;
    //cout << "potential energy "<< potentialenergy<<endl;

    for(unsigned int i=0; i<crystall->allatoms.size(); i++){
        updateVelocity(crystall->allatoms[i]);

        //kinetic energy of the atom in the crystal
        crystall->energy+=0.5*dot(crystall->allatoms[i]->getVelocity(),crystall->allatoms[i]->getVelocity());
    }
    crystall->ke= crystall->energy-crystall->pe;
// cout <<"crystal energy "<<crystall->energy <<endl;
    if(thermalize){
        thermostatBerendsen();
    }
    if(crystall->beginenergy==0&&!thermalize){
        crystall->beginenergy=crystall->energy;
    }
}

void VerletAlgo2::thermostatAnders(){
    double tem=crystall->temperature();
    double ratio = crystall->inittemp/tem;
    for(int i=0; i<crystall->allatoms.size();i++){
        Atom *atom = crystall->allatoms[i];
        vec3 velocity = atom->getVelocity();
        atom->setVelocity(velocity*ratio);
    }
}

void VerletAlgo2::thermostatBerendsen(){
    double tau=20.0*h;
    double tem=crystall->temperature();
    double tbath=crystall->inittemp;
    double gamma = 1+h/tau*(tbath/tem-1);
    gamma=sqrt(gamma);
    for(int i=0; i<crystall->allatoms.size();i++){
        Atom *atom = crystall->allatoms[i];
        vec3 velocity = atom->getVelocity();
        atom->setVelocity(velocity*gamma);
    }
}

void VerletAlgo2::thermostatAndersen(){
    double tau=20.0*h;
    int nAtoms= crystall->numberofatoms;
    vec  v = randu<vec>(nAtoms);
    double tem=this->crystall->inittemp;
    for(int i=0; i<nAtoms;i++){
        if(v(i)<h/tau){
            Atom *atom = crystall->allatoms[i];
            vec3 v1; v1.zeros();
            v1 << DRanNormalZigVec()*sqrt(tem)<< DRanNormalZigVec()*sqrt(tem)<< DRanNormalZigVec()*sqrt(tem);
            atom->setVelocity(v1);
        }
    }
//
}

void VerletAlgo2::integrate_noapprox(){
    //ofstream debugging;
    //debugging.open("/home/jonathan/projectsFSAP/project1/project1/debuglog2.txt", ios::app);
    for(unsigned int i=0; i<crystall->allatoms.size(); i++){
        updateVelocity(crystall->allatoms[i]);
        updatePosition(crystall->allatoms[i]);
        vec3 nullvector; nullvector.zeros();
        crystall->allatoms[i]->setAcceler(nullvector);
    }
// cout << crystall->temperature()<<endl;
    for(unsigned int i=0; i<crystall->allatoms.size(); i++){
        updateAccelerNoApprox(crystall->allatoms[i]);

    }
    for(unsigned int i=0; i<crystall->allatoms.size(); i++){
        updateVelocity(crystall->allatoms[i]);
    }
}

void VerletAlgo2::updateAcceler(Atom *atom){
    bool debugg=false;

    Atom* otheratom=atom->nextAtom;
    while(otheratom!=NULL){
        calcForce(atom, otheratom);
        otheratom=otheratom->nextAtom;
    }

    //ofstream debugging;
    //debugging.open("/home/jonathan/projectsFSAP/project1/project1/debuglog2.txt", ios::app);
    //indices of cell atom is in nrXYZ:
    int nrXYZ[3];
    vec3 r= atom->getPosition();
    for(int i=0; i<3; i++){
        nrXYZ[i]=int(r(i)/crystall->vectorBC(i));
    }
    int nrX[3], nrY[3], nrZ[3];
    findXYZCellIndices(nrXYZ, nrX, nrY, nrZ);
    //indices of all neighbouring cells are now in nrX, nrY and nrZ

    for(int i=0; i<3; i++){
        for(int j=0; j<3;j++){
            for(int k=0;k<3;k++){
                if(nrX[i]!=-1 &nrY[j]!=-1&nrZ[k]!=-1){
                    if(crystall->allcells.at(nrX[i]).at(nrY[j]).at(nrZ[k]).visited!=true){
                        otheratom = crystall->allcells.at(nrX[i]).at(nrY[j]).at(nrZ[k]).first;
                        while(otheratom!=NULL){
                            calcForce(atom, otheratom);
                            otheratom=otheratom->nextAtom;
                        }
                    }
                }
            }
        }
    }
}

void VerletAlgo2::updateAccelerNoApprox(Atom *atom){
    bool debugg=false;
    //ofstream debugging;
    //debugging.open("/home/jonathan/projectsFSAP/project1/project1/debuglog2.txt", ios::app);
    //indices of cell atom is in nrXYZ:
// int nrXYZ[3];
// vec3 r= atom->getPosition();
// for(int i=0; i<3; i++){
// nrXYZ[i]=int(r(i)/crystall->vectorBC(i));
// }
// int nrX[3], nrY[3], nrZ[3];
// findXYZCellIndices(nrXYZ, nrX, nrY, nrZ);
    //indices of all neighbouring cells are now in nrX, nrY and nrZ

    for(int i=atom->number; i<crystall->allatoms.size();i++){
        Atom* otheratom = crystall->allatoms[i];
        calcForce(atom, otheratom);
    }

}

void VerletAlgo2::calcForce(Atom* atom, Atom* otheratom){
    int i = atom->number;
    int j = otheratom->number;

    //ofstream debugging;
    //debugging.open("/home/jonathan/projectsFSAP/project1/project1/debuglog2.txt", ios::app);
    //stop if integrating atom with itself, should normally not happen
    if(i==j){
        return;
    }

    vec3 position=atom->getPosition();
    vec3 othervec=otheratom->getPosition();
    vec3 closestvector = findClosestPosition(position, othervec);
    vec3 relvec = position - closestvector;
    vec3 relvec2 = position - othervec;
    double r2=dot(relvec,relvec);
    double r6=r2*r2*r2;
    double r12=r6*r6;

    vec3 oneacceler=atom->getAcceler();
    vec3 otheracceler=otheratom->getAcceler();
    for(int k=0; k<3; k++){
        double temp = 24*(2.0/r12-1.0/r6)*relvec(k)/r2;
        crystall->energy+=2.0*LJpotential( relvec2);
        if(temp>cutoffacceleration){
            temp=cutoffacceleration;
            //crystall->forces(i,j,k)=cutoffacceleration;
            //crystall->forces(j,i,k)=-1*cutoffacceleration;
        }
        else if(temp<-1*cutoffacceleration){
            temp=-cutoffacceleration;
        }
// if(atom->number==0 && otheratom->number==586){
// cout << "posatom0" << position <<endl;
// cout << "posatom586"<<othervec<<endl;
// cout << "closestvector"<<closestvector<<endl;
// cout << relvec << endl;
// cout << temp <<endl;
// cout << "r2 "<<r2 << " r6 " <<r6 << " r12 " << r12 <<endl;
// cout <<24*(2.0/r12-1.0/r6)/r2<<endl;
// }

        oneacceler(k)+=temp;
        otheracceler(k)-=temp;

        crystall->pressure+=temp*relvec(k);
    }
    atom->setAcceler(oneacceler);
    otheratom->setAcceler(otheracceler);
}

//this function finds the cellindices of the neighbouring cells
//it takes the minimal image convention into count
//it puts all X, Y and Z indices of the neighbouring cells in the last 3
//function arguments
//!!!!!!!!!!!!!!!!!function is debugged!!!!!!!!!!!!
void VerletAlgo2::findXYZCellIndices(int* nrXYZ, int* nrX, int* nrY, int* nrZ){
    //the maximum indices of cells in x, y and z direction
    int imax = crystall->allcells.size();
    int jmax=crystall->allcells.at(0).size();
    int kmax=crystall->allcells.at(0).at(0).size();

    int i=0;
    for(int l=nrXYZ[0]-1; l<nrXYZ[0]+2; l++){
        if(l<0){
            nrX[i]=l+imax;
        }
        else if(l>=imax){
            nrX[i]=l-imax;
        }
        else{
            nrX[i]=l;
        }
        i++;
    }
    i=0;
    for(int l=nrXYZ[1]-1; l<nrXYZ[1]+2; l++){
        if(l<0){
            nrY[i]=l+jmax;
        }
        else if(l>=jmax){
            nrY[i]=l-jmax;
        }
        else{
            nrY[i]=l;
        }
        i++;
    }
    i=0;
    for(int l=nrXYZ[2]-1; l<nrXYZ[2]+2; l++){
        if(l<0){
            nrZ[i]=l+kmax;
        }
        else if(l>=kmax){
            nrZ[i]=l-kmax;
        }
        else{
            nrZ[i]=l;
        }
        i++;
    }
    for(i=0; i<3;i++){
        for(int j=i+1;j<3; j++){
            if(nrX[i]==nrX[j]){
// cout<< "lalala" <<endl;
// cout<< "i "<<i << " j "<<j << " nrX "<< nrX[i]<<endl;
                nrX[j]=-1;
            }
            if(nrY[i]==nrY[j]){
// cout<< "lalala" <<endl;
// cout<< "i "<<i << " j "<<j << " nrY "<< nrY[i]<<endl;
                nrY[j]=-1;
            }
            if(nrZ[i]==nrZ[j]){
// cout<< "lalala" <<endl;
// cout<< "i "<<i << " j "<<j << " nrZ "<< nrZ[i]<<endl;
                nrZ[j]=-1;
            }
        }
    }

}

vec3 VerletAlgo2::findClosestPosition(vec3 &position, vec3 &otherposition){
    vec3 answer; answer.fill(0);
    bool debugg=false;
    for(int i=0; i<3; i++){
        double projectionother = otherposition(i);
        double projectionpos = position(i);
        double l=this->crystall->boundary(i);
        double distance=l;
        for(int j=-1; j<2; j++){
            distance=abs(projectionpos-(projectionother+j*l));
            if(distance<=l/2){
                answer(i)=projectionother+j*l;
            }
        }
    }
    return answer;
}

void VerletAlgo2::updateVelocity(Atom *atom){
    vec3 velocity=atom->getVelocity();
    vec3 acceler=atom->getAcceler();
    velocity+=0.5*acceler*this->h;
    atom->setVelocity(velocity);
}

void VerletAlgo2::updatePosition(Atom *atom){
    vec3 position=atom->getPosition();
    vec3 velocity=atom->getVelocity();

    int nrXYZ[3];
    for(int i=0; i<3; i++){
        nrXYZ[i]=int(position(i)/crystall->vectorBC(i));
    }
    position+=velocity*this->h;

    atom->realposition+=velocity*this->h;
    vec3 finalposition =boundCheck(position);
    atom->setPosition(finalposition);
    for(int i=0; i<3; i++){
        this->crystall->msqdplm+=(atom->realposition(i)-atom->initialposition(i))*(atom->realposition(i)-atom->initialposition(i));
    }

    //if atom is not anymore in its cell
    if(!atom->currentcell->isAtomInCell(atom)){
        atom->currentcell->removeelement(atom);
        int x,y,z;
        crystall->findCellOfAtom(atom, x, y, z);
        crystall->allcells.at(x).at(y).at(z).insertElement(atom);
    }
}


vec3 VerletAlgo2::boundCheck(vec3 &position){
    vec3 boundvec = this->crystall->boundary;
    vec3 answer=position;
    for(int i=0; i<3; i++){
        while(answer(i)<0){
            //debugging << "summing boundvec at i="<<i<< " position(i)=" << position(i);
            answer(i)+=boundvec(i);
            //debugging << "new position(i) " << position(i) << endl;
        }
        while(answer(i)>=boundvec(i)){
            answer(i)-=boundvec(i);
        }
    }
    return answer;
}

double VerletAlgo2::LJpotential(vec3 relvec){
    double r2=relvec(0)*relvec(0)+relvec(1)*relvec(1)+relvec(2)*relvec(2);
    double r6=r2*r2*r2;
    double r12=r6*r6;

    return 4*(1.0/r12-1.0/r6);
}
