#include <iostream>
#include <fstream>
#include "crystal.h"
#include <armadillo>
#include <ctime>
#include <vector>
#include "verletalgo2.h"
#include <sstream>
#include "printing.h"
#include <mpi.h>


using namespace arma;
using namespace std;


void calculatepressures(double _temp, double *_avgpres, double *_stddevpres, double _b){
    time_t tinit = time(0);



    int seed = -1;

//Liquid:
//T=120K, rho/rho_0=0.8
//b=5.82A
//P=38MPa

//SOLID:
//T=60K, rho/rho_0=1.2
//b=5.09A
//P=0.6GPa

//Gas:
//T=359K
//rho/rho_0 = 0.3
//b=8.07A
//P=44MPa


    int nc=8;
    double h=0.01;
    //latice parameter in unit Angstrom
    double b=_b;
    double temperature=_temp;

    //200
    int nrofthermalizingsteps=200;

    int nrofstepstotakemeasurements=1000;

    system("rm /home/jonathan/projectsFSAP/project1/project1/output/*.xyz");
    Printing p;
    Crystal crystal(nc, b, seed, temperature);
    p.printing(crystal);
    //VerletAlgo integrator(crystal);
    VerletAlgo2 integrator(&crystal, h);

//    ofstream measurements;
//    measurements.open("/home/jonathan/projectsFSAP/project1/project1/output/measurements-gas.txt");


    double pressures[nrofstepstotakemeasurements];

//    measurements << "#1:timesteps #2: time #3:temperature  #4:total energy #5:kinetic energy #6:potential energy #7:relative energy error #8:pressure"<<endl;
    for(int j=0; j<nrofthermalizingsteps; j++){
        integrator.integrate(true);
//        measurements << j<<" "<< j*integrator.h << " " <<crystal.temperature()<<" " << integrator.crystall->energy<<" " << integrator.crystall->ke<<" " << integrator.crystall->pe << " " << abs((integrator.crystall->energy - integrator.crystall->beginenergy)/integrator.crystall->beginenergy) << endl;
        if(j%100==0){
            cout << "now in step " << j << " in the thermilisation phase" << endl;
            //cout << "nrofatomsfound "<<crystal.countAtoms()<< endl;
        }
        //ofstream output;
        //output.open(p.createname(j).c_str());

        //output << crystal << endl;
    }
    for(int j=0; j<nrofstepstotakemeasurements; j++){

        integrator.integrate(false);


        //integrator.integrate_noapprox();
        pressures[j]=integrator.crystall->pressure;
//        measurements << j<<" "<< j*integrator.h << " " <<crystal.temperature()<<" " << integrator.crystall->energy<<" " << integrator.crystall->ke<<" " << integrator.crystall->pe << " " << abs((integrator.crystall->energy - integrator.crystall->beginenergy)/integrator.crystall->beginenergy) <<" "<< integrator.crystall->pressure<< endl;
        //p.printvelocities(crystal, j);
        if(j%50==0){
//            cout << "now in step " << j << " doing measurements, i am processor  " << myrank << endl;
        }
    }
    //int j=0;
    //p.printvelocities(crystal, j);
    double avgpres=0;
    for(int j=0; j<nrofstepstotakemeasurements; j++){
        avgpres+=pressures[j]/nrofstepstotakemeasurements;
    }

    *_avgpres=avgpres;

//    cout << "average final temperature " << avgtemp << endl;

    double stddev=0;
    for(int j=0; j<nrofstepstotakemeasurements; j++){
        stddev+=(1.0/(nrofstepstotakemeasurements-1))*(pressures[j]-avgpres)*(pressures[j]-avgpres);
    }
    stddev=sqrt(stddev);

    *_stddevpres=stddev;


//    cout << "stddev final temperature " << stddev << endl;

//    cout << "the temperature is " << integrator.crystall->temperature() << "for processor " << myrank <<endl;
//    cout << "the boundary vector is " << crystal.boundary << endl;
//    cout << "BC vector is " << crystal.vectorBC << endl;
//    cout << "crystal has " << integrator.crystall->countAtoms() << " atoms "<< endl;
    time_t tdone = time(0);
//    cout << "I AM PROCESSOR "<< myrank << " and i worked for " << tdone-tinit << " seconds " << endl;


}

string makeDirectory(){
    stringstream oss;
    oss<<"mkdir /home/jonathan/projectsFSAP/project2/output/";
    time_t t= time(0);
    oss<<t;
    system(oss.str().c_str());
    stringstream answer;
    answer << "/home/jonathan/projectsFSAP/project2/output/";
    answer << t<<"/";
    return answer.str();
}

string outputFile(const string dirname, const string filename){
    stringstream answer;
    answer<< dirname<<filename;
    return answer.str();
}

void simulationwithoutput(){
    Printing p;

    time_t tinit = time(0);
    string dirname = makeDirectory();

    //Liquid:
    //T=120K, rho/rho_0=0.8
    //b=5.82A
    //P=38MPa

    int seed = -1;
    int nc=20;
    double h=0.005;
    //latice parameter in unit Angstrom
    double b=5.72;
    double temperature=0.851*tempunit;
    int nrofthermalizingsteps=500;
    int nrofstepstotakemeasurements=0;



//    system("rm /home/jonathan/projectsFSAP/project1/project1/output/*.xyz");
//    Crystal crystal(nc, b, seed, temperature);

    string name="/home/jonathan/projectsFSAP/project2/output/1363095443/locationatoms.00499.xyz";
    Crystal crystal(name);
    for(int i=0; i<20; i++){
        crystal.createSphere(20.0, 30.0);
        string filename = p.createname(i);
        string completename = outputFile(dirname, filename);
        ofstream output;
        output.open(completename.c_str());
        output << crystal << endl;
    }



    return;

    //VerletAlgo integrator(crystal);
    VerletAlgo2 integrator(&crystal, h);





    ofstream measurements;
    string measname = outputFile(dirname, "measurements.txt");
    measurements.open(measname.c_str());
    measurements << "#1:timesteps #2: time #3:temperature  #4:total energy #5:kinetic energy #6:potential energy #7:relative energy error #8:pressure #9:mean square displacement"<<endl;
    for(int j=0; j<nrofthermalizingsteps; j++){
        integrator.integrate(true);
        string filename = p.createname(j);
        string completename = outputFile(dirname, filename);
        measurements << j<<" "<< j*integrator.h << " " <<crystal.temperature()<<" " << integrator.crystall->energy<<" " << integrator.crystall->ke<<" " << integrator.crystall->pe << " " << abs((integrator.crystall->energy - integrator.crystall->beginenergy)/integrator.crystall->beginenergy) <<" "<< integrator.crystall->pressure<< " " << integrator.crystall->msqdplm<<endl;

        if(j%50==0){
            cout << "now in step " << j << " in the thermilisation phase" << endl;
        }
        ofstream output;
        output.open(completename.c_str());
        output << crystal << endl;

    }
    return;

    for(int j=0; j<nrofstepstotakemeasurements; j++){

        integrator.integrate(false);


        //integrator.integrate_noapprox();
        measurements << j<<" "<< j*integrator.h << " " <<crystal.temperature()<<" " << integrator.crystall->energy<<" " << integrator.crystall->ke<<" " << integrator.crystall->pe << " " << abs((integrator.crystall->energy - integrator.crystall->beginenergy)/integrator.crystall->beginenergy) <<" "<< integrator.crystall->pressure<< " " << integrator.crystall->msqdplm<<endl;
        if(j%50==0){
            cout << "now in step " << j << " doing measurements, i am processor  "  << endl;
        }
        ofstream output;
//        output.open(p.createname(j).c_str());

        output << crystal << endl;

    }
//    integrator.crystall->radialDistFunction();


}

void normalsimulation(){
    time_t tinit = time(0);

    //Liquid:
    //T=120K, rho/rho_0=0.8
    //b=5.82A
    //P=38MPa

    int seed = -1;
    int nc=7;
    double h=0.0025;
    //latice parameter in unit Angstrom
    double b=5.82;
    double temperature=tempunit;

    //200
    int nrofthermalizingsteps=200;

    int nrofstepstotakemeasurements=1000;

    system("rm /home/jonathan/projectsFSAP/project1/project1/output/*.xyz");
    Crystal crystal(nc, b, seed, temperature);

    //VerletAlgo integrator(crystal);
    VerletAlgo2 integrator(&crystal, h);

    for(int j=0; j<nrofthermalizingsteps; j++){
        integrator.integrate(true);
        if(j%50==0){
            cout << "now in step " << j << " in the thermilisation phase" << endl;
        }
    }
    for(int j=0; j<nrofstepstotakemeasurements; j++){

        integrator.integrate(false);
        if(j%50==0){
            cout << "now in step " << j << " taking measurement" << endl;
        }
    }


//    cout << "stddev final temperature " << stddev << endl;

//    cout << "the temperature is " << integrator.crystall->temperature() << "for processor " << myrank <<endl;
//    cout << "the boundary vector is " << crystal.boundary << endl;
//    cout << "BC vector is " << crystal.vectorBC << endl;
//    cout << "crystal has " << integrator.crystall->countAtoms() << " atoms "<< endl;
    time_t tdone = time(0);
    cout << "I AM PROCESSOR " << " and i worked for " << tdone-tinit << " seconds " << endl;


}

void pressuresimulation(){
    //UNIT SYSTEM!!
    //Distances = Angstrom
    //Time = picosecond
    //Energy = electronvolt
    //Mass = enter in amu, wrong unit system is addressed for in value of k
    //Temperature = Kelvin
    //velocity = angstrom/ps
    //k = 0.8314766196505026 when masses are expressed in amu and T in K to get velocities in A/ps

    //cubic lattice with nc x nc x nc cells

    int argc; char** argv; int numprocs, myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    ofstream presstream;
    presstream.open("/home/jonathan/projectsFSAP/project1/project1/output/pressuresdensity.txt", ios::app);

    double pres, presstdv;
    for(int j=3; j<10; j++){
        for(int i=0; i<40; i++){
            if(myrank==i%numprocs){
                double b=double(j);
                double density = 4.0*8.0*8.0*8.0/(xunit*xunit*xunit);
                double temperature=i*20;
                cout << "i am processor "<< myrank << " and i start simulation with temperature T="<<temperature << "and b=" << b <<endl;
                calculatepressures(temperature, &pres, &presstdv, b);
                presstream << temperature/tempunit << " " << density << pres << " " << presstdv <<endl;
                cout << "i am processor " << myrank<< " " <<temperature << " " << pres << " " << presstdv <<endl;
            }
        }

    }




    MPI_Finalize();
}

int main()
{
    simulationwithoutput();
    return 0;
}


