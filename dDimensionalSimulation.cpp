#include "std_include.h" // std library includes, definition of scalar, etc.. has a "using namespace std" in it, because I'm lazy

//we'll use TCLAP as our command line parser
#include <tclap/CmdLine.h>
#include "cuda_profiler_api.h"

#include "functions.h"
#include "gpuarray.h"
#include "periodicBoundaryConditions.h"
#include "simulation.h"
#include "simpleModel.h"
#include "baseUpdater.h"
#include "energyMinimizerFIRE.h"
#include "velocityVerlet.h"
#include "noseHooverNVT.h"
#include "noiseSource.h"
#include "harmonicRepulsion.h"
#include "lennardJones6_12.h"
#include "indexer.h"
#include "hyperrectangularCellList.h"
#include "neighborList.h"
#include "kdTreeNeighborList.h"
#include "poissonDiskSampling.h"

using namespace std;
using namespace TCLAP;

//!What, after all, *is* the volume of a d-dimensional sphere?
scalar sphereVolume(scalar radius, int dimension)
    {
    if(dimension == 1)
        return 2*radius;
    else
        if(dimension == 2)
            return PI*radius*radius;
        else
            return (2.*PI*radius*radius)/((scalar) dimension)*sphereVolume(radius,dimension-2);
    };

double estimateEnergy(shared_ptr<simpleModel> conf, shared_ptr<kdTreeNeighborList> nl, shared_ptr<force> f)
    {
    shared_ptr<Simulation> sim = make_shared<Simulation>();
    sim->setConfiguration(conf);
    sim->setBox(conf->Box);
    f->setNeighborList(nl);
    sim->addForce(f,conf);
    return sim->computePotentialEnergy();
    }

void estimateTrueEnergy(shared_ptr<simpleModel> conf, shared_ptr<baseNeighborList> nl, shared_ptr<kdTreeNeighborList> nlEst,shared_ptr<force> f,double &e, double &eEst)
    {
    shared_ptr<Simulation> sim = make_shared<Simulation>();
    sim->setConfiguration(conf);
    sim->setBox(conf->Box);
    f->setNeighborList(nlEst);
    sim->addForce(f,conf);
    eEst =  sim->computePotentialEnergy();

    sim->clearForceComputers();
    f->setNeighborList(nl);
    sim->addForce(f,conf);
    e =  sim->computePotentialEnergy();
    }

/*!
This file runs some dynamics on particles interacting according to some
potential... when this repository is ever meant to be used this all should be
updated.
*/
int main(int argc, char*argv[])
{
    // wrap tclap in a try block
    try
    {
    //First, we set up a basic command line parser...
    //      cmd("command description message", delimiter, version string)
    CmdLine cmd("basic testing of dDimSim", ' ', "V0.1");

    //define the various command line strings that can be passed in...
    //ValueArg<T> variableName("shortflag","longFlag","description",required or not, default value,"value type",CmdLine object to add to
    ValueArg<int> programSwitchArg("z","programSwitch","an integer controlling program branch",false,0,"int",cmd);
    ValueArg<int> gpuSwitchArg("g","USEGPU","an integer controlling which gpu to use... g < 0 uses the cpu",false,-1,"int",cmd);
    ValueArg<int> nSwitchArg("n","Number","number of particles in the simulation",false,100,"int",cmd);
    ValueArg<int> initialIterationsSwitchArg("i","iterations","number of timestep iterations",false,100,"int",cmd);
    ValueArg<int> maxIterationsSwitchArg("j","otherIterations","number of timestep iterations",false,100,"int",cmd);
    ValueArg<scalar> lengthSwitchArg("l","sideLength","size of simulation domain",false,10.0,"double",cmd);
    ValueArg<scalar> temperatureSwitchArg("t","temperature","temperature of simulation",false,.00,"double",cmd);

    //allow setting of system size by either volume fraction or density (assuming N has been set)
    scalar phiDest = 1.90225*exp(-(scalar)DIMENSION / 2.51907);
    ValueArg<scalar> phiSwitchArg("p","phi","volume fraction",false,phiDest,"double",cmd);
    ValueArg<scalar> rhoSwitchArg("r","rho","density",false,-1.0,"double",cmd);
    //parse the arguments
    cmd.parse( argc, argv );

    int programSwitch = programSwitchArg.getValue();
    int N = nSwitchArg.getValue();
    int maximumIterations = initialIterationsSwitchArg.getValue();
    int iterations = maxIterationsSwitchArg.getValue();
    scalar L = lengthSwitchArg.getValue();
    scalar Temperature = temperatureSwitchArg.getValue();
    scalar phi = phiSwitchArg.getValue();
    scalar rho = rhoSwitchArg.getValue();

    int gpuSwitch = gpuSwitchArg.getValue();
    bool GPU = false;
    if(gpuSwitch >=0)
        GPU = chooseGPU(gpuSwitch);

    if(phi >0)
        {
        L = pow(N*sphereVolume(.5,DIMENSION) / phi,(1.0/(scalar) DIMENSION));
        rho = N/pow(L,(scalar)DIMENSION);
        }
    else
        phi = N*sphereVolume(.5,DIMENSION) / pow(L,(scalar)DIMENSION);

    if(rho >0)
        {
        L = pow(((scalar)N/rho),(1.0/(scalar) DIMENSION));
        phi = rho * sphereVolume(.5,DIMENSION);
        }
    else
        rho = N/pow(L,(scalar)DIMENSION);


    int dim =DIMENSION;
    cout << "running a simulation in "<<dim << " dimensions with box sizes " << L << endl;
    cout << "density = " << rho << "\tvolume fracation = "<<phi<<endl;
    shared_ptr<simpleModel> Configuration = make_shared<simpleModel>(N);
    shared_ptr<periodicBoundaryConditions> PBC = make_shared<periodicBoundaryConditions>(L);

    shared_ptr<Simulation> sim = make_shared<Simulation>();
    sim->setConfiguration(Configuration);
    sim->setBox(PBC);

    noiseSource noise(true);
    Configuration->setParticlePositionsRandomly(noise);
    scalar ke = Configuration->setVelocitiesMaxwellBoltzmann(Temperature,noise);
    printf("temperature input %f \t temperature calculated %f\n",Temperature,Configuration->computeInstantaneousTemperature());

    double range = 1.0;
    double epsilon = 0.0;
    shared_ptr<neighborList> neighList = make_shared<neighborList>(range,PBC,1);
    shared_ptr<kdTreeNeighborList> kdNeighList = make_shared<kdTreeNeighborList>(range,PBC,epsilon,true);


     //monodisperse harmonic spheres
    shared_ptr<harmonicRepulsion> softSpheres = make_shared<harmonicRepulsion>();
    softSpheres->setMonodisperse();
    softSpheres->setNeighborList(neighList);
    vector<scalar> stiffnessParameters(1,1.0);
    softSpheres->setForceParameters(stiffnessParameters);
    sim->addForce(softSpheres,Configuration);
    cout << "simulation set-up finished" << endl;cout.flush();

    shared_ptr<noseHooverNVT> nvt = make_shared<noseHooverNVT>(Configuration,Temperature);
    nvt->setDeltaT(1e-2);
    //shared_ptr<energyMinimizerFIRE> fire = make_shared<energyMinimizerFIRE>(Configuration);
    //fire->initializeParameters();
    //fire->setForceCutoff(1e-10);
    sim->addUpdater(nvt,Configuration);

    if(gpuSwitch >=0)
        {
        sim->setCPUOperation(false);
//        Configuration->setGPU();
//        softSpheres->setGPU();
//        fire->setGPU();
//        neighList->setGPU();
        };

    clock_t t1 = clock();
    cudaProfilerStart();

    scalar dt=-12;
    for (int timestep = 0; timestep < maximumIterations; ++timestep)
        {
        sim->performTimestep();
        if(timestep%100 == 0)
            printf("timestep %i: target T = %f\t instantaneous T = %g\t PE = %g\t nlist max = %i\n",timestep,Temperature,Configuration->computeInstantaneousTemperature(),sim->computePotentialEnergy(),neighList->Nmax);
        };
    cudaProfilerStop();
    clock_t t2 = clock();

    cout << endl << endl <<  "intialization done" << endl << endl;

    shared_ptr<kdTreeNeighborList> kdNeighs = make_shared<kdTreeNeighborList>(range,PBC,epsilon,true);

    double pe;
    double peEst;
    double peLast=0;
    double peLastEst=0;
    double rangeEst = 1.;
    double epsilonEst = .1;

    char filename[256];
    sprintf(filename,"../data/test_N%i.txt",N);
    ofstream myfile;
    myfile.open(filename);

    sim->setCPUOperation(true);

    for (int timestep = 0; timestep < iterations; ++timestep)
        {
        if(timestep%100==0)
            {
            printf("timestep %i: target T = %f\t instantaneous T = %g\t PE = %g\t nlist max = %i\n",timestep,Temperature,Configuration->computeInstantaneousTemperature(),sim->computePotentialEnergy(),neighList->Nmax);
            }
        sim->performTimestep();
        //printf("timestep %i: target T = %f\t instantaneous T = %g\t PE = %g\t nlist max = %i\n",timestep,Temperature,Configuration->computeInstantaneousTemperature(),sim->computePotentialEnergy(),neighList->Nmax);
        shared_ptr<kdTreeNeighborList> kdNeighEst = make_shared<kdTreeNeighborList>(rangeEst,PBC,epsilonEst,true);
        estimateTrueEnergy(Configuration,neighList,kdNeighEst,softSpheres,pe,peEst);

        cout << "true energy = " << pe << " estimated = " << peEst << endl;
        cout<< " DeltaE - DeltaEEst = " << (pe-peLast)- (peEst-peLastEst) << endl;
        if(timestep > 0)
            myfile << (pe-peLast)- (peEst-peLastEst) << "\t" << pe << "\t" << peLast <<"\n";
        peLast = pe;
        peLastEst=peEst;
        };

    myfile.close();

    scalar E = sim->computePotentialEnergy();
    printf("simulation potential energy at %f\n",E);

    scalar timeTaken = (t2-t1)/(scalar)CLOCKS_PER_SEC/maximumIterations;
    cout << endl << "simulations took " << timeTaken << " per time step" << endl << endl;

/*
    ofstream ofs;
    char dataname[256];
    sprintf(dataname,"../data/timing_d%i_g%i.txt",DIMENSION,gpuSwitch);
    ofs.open(dataname,ofstream::app);
    ofs << N <<"\t" << timeTaken << "\n";
    ofs.close();
    t1 = clock();
    neighList->computeNeighborLists(Configuration->returnPositions());
    t2 = clock();
    scalar ntime = (t2-t1)/(scalar)CLOCKS_PER_SEC;
    cout << endl << "nlists take " << ntime << endl;
    t1 = clock();
    softSpheres->computeForces(Configuration->returnForces());
    t2 = clock();
    scalar ftime = (t2-t1)/(scalar)CLOCKS_PER_SEC - ntime;
    cout << endl << "forces take " << ftime << endl;
    t1 = clock();
    nve->performUpdate();
    t2 = clock();
    scalar stime = (t2-t1)/(scalar)CLOCKS_PER_SEC - ntime - ftime;
    cout << endl << "timestep takes" << stime << endl;

*/

//
//The end of the tclap try
//
    } catch (ArgException &e)  // catch any exceptions
    { cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }
    return 0;
};
