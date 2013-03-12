#include "printing.h"

Printing::Printing()
{
}

void Printing::printing(Crystal &crystal){
    ofstream output;
    output.open("/home/jonathan/projectsFSAP/project2/project2/output/locationatoms.xyz");
    output << crystal << endl;
    output.close();
}

void Printing::printvelocities(Crystal &crystal, int j){
    int number = j;
    stringstream oss;
    oss<<"/home/jonathan/projectsFSAP/project2/project2/output/";
    oss<<"velocity-vtot-.";
    if(number<10){
        oss<<"000"<<number;
    }
    else if(number <100){
        oss<<"00"<<number;
    }
    else if(number <1000){
        oss<<"0"<<number;
    }
    else{
        oss<<number;
    }

    oss<<".dat";

    ofstream outputx, outputy, outputz, outputtot;
    //outputx.open(();
    //outputx.open(oss.str().c_str());
    //outputy.open("/home/jonathan/projectsFSAP/project2/project2/output/velocity-y.dat");
    //outputz.open("/home/jonathan/projectsFSAP/project2/project2/output/velocity-z.dat");
    outputtot.open(oss.str().c_str());

    for(unsigned int i=0; i<crystal.allatoms.size(); i++){
        vec3 velocity= crystal.allatoms[i]->getVelocity();
        outputtot << norm(crystal.allatoms[i]->getVelocity(),2) << endl;
        //outputx << velocity(0) << endl;
        //outputy << velocity(1) << endl;
        //outputz << velocity(2)<< endl;
    }
    outputx.close();
    //outputy.close();
    //outputz.close();
    outputtot.close();
}

string Printing::createname( int number){
    stringstream oss;
//    oss<<"/home/jonathan/projectsFSAP/project2/project2/output/";
    oss<<"locationatoms.";
    oss<<setfill('0')<<setw(5)<<number<<".xyz";
    return oss.str();
    if(number<10){
        oss<<"0000"<<number;
    }
    else if(number <100){
        oss<<"000"<<number;
    }
    else if(number <1000){
        oss<<"00"<<number;
    }
    else if(number <10000){
        oss <<"0"<<number;
    }
    else{
        oss<<number;
    }
    oss<<".xyz";
//    return oss.str();
}
