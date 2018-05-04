#include "std_include.h" // std library includes, definition of scalar, etc.. has a "using namespace std" in it, because I'm lazy

//we'll use TCLAP as our command line parser
#include <tclap/CmdLine.h>
#include "functions.h"
#include "gpuarray.h"
#include "periodicBoundaryConditions.h"
#include "simulation.h"
#include "simpleModel.h"
#include "baseUpdater.h"
#include "energyMinimizerFIRE.h"
#include "noiseSource.h"
#include "harmonicBond.h"
#include "harmonicAngle.h"

#include "indexer.h"
#include "hyperrectangularCellList.h"

using namespace std;
using namespace TCLAP;

/*!
command line parameters help identify a data directory and a filename... the output is a text file
(in the data/ directory rooted here) containing easy-to-read fourier transforms of the height-map
representation of the extremal interfaces for each point in time
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
    ValueArg<scalar> lengthSwitchArg("l","sideLength","size of simulation domain",false,10.0,"double",cmd);
    //parse the arguments
    cmd.parse( argc, argv );

    int programSwitch = programSwitchArg.getValue();
    int N = nSwitchArg.getValue();
    scalar L = lengthSwitchArg.getValue();
    int gpuSwitch = gpuSwitchArg.getValue();
    bool GPU = false;
    if(gpuSwitch >=0)
        GPU = chooseGPU(gpuSwitch);

    int dim =DIMENSION;
    cout << "running a simulation in "<<dim << " dimensions" << endl;
    shared_ptr<simpleModel> Configuration = make_shared<simpleModel>(N);
    shared_ptr<periodicBoundaryConditions> PBC = make_shared<periodicBoundaryConditions>(L);
    shared_ptr<updater> upd = make_shared<updater>(1);
    shared_ptr<energyMinimizerFIRE> fire = make_shared<energyMinimizerFIRE>(Configuration);
    fire->setFIREParameters(0.02,0.99,0.1,1.1,0.95,.9,4,1e-12);
    fire->setMaximumIterations(1000);

    shared_ptr<hyperrectangularCellList> cellList = make_shared<hyperrectangularCellList>(1.5,PBC);

    shared_ptr<harmonicBond> bonds = make_shared<harmonicBond>();
    vector<simpleBond> blist;
    simpleBond testBond(0,1,1.0,1.0);
    blist.push_back(testBond);
    simpleBond testBond2(2,1,4.,1.0);
    blist.push_back(testBond2);
    testBond.setBondIndices(5,6);
    testBond.setRestLength(1.02);
    blist.push_back(testBond);
    bonds->setBondList(blist);

    shared_ptr<harmonicAngle> angles = make_shared<harmonicAngle>();
    vector<simpleAngle> alist;
    simpleAngle testAngle(0,1,2,PI/3.0,1.0);
    alist.push_back(testAngle);
    angles->setAngleList(alist);


    shared_ptr<Simulation> sim = make_shared<Simulation>();
    sim->setConfiguration(Configuration);
    sim->setBox(PBC);

    //after the simulation box has been set, we can set particle positions
    noiseSource noise(true);
    Configuration->setParticlePositionsRandomly(noise);

    /*
    //for testing, print positions
    {
    ArrayHandle<dVec> pos(Configuration->returnPositions());
    for (int pp = 0; pp < N; ++pp)
        printdVec(pos.data[pp]);
    }
    */

    cout << endl << endl;
    sim->addForce(bonds,Configuration);
    sim->addForce(angles,Configuration);
    sim->addUpdater(fire,Configuration);
    clock_t t1 = clock();
    sim->performTimestep();
    clock_t t2 = clock();

    /*
    //how did FIRE do? check by hand
    {
    ArrayHandle<dVec> pos(Configuration->returnPositions());
    for (int pp = 0; pp < N; ++pp)
        printdVec(pos.data[pp]);
    }
    */
    cout << endl << "minimization took " << (t2-t1)/(scalar)CLOCKS_PER_SEC << endl;
    if(gpuSwitch >=0) cellList->setGPU(true);
    t1 = clock();
    cellList->computeCellList(Configuration->returnPositions());
    t2 = clock();
    cout << endl << "cellList took " << (t2-t1)/(scalar)CLOCKS_PER_SEC << " and iterated through the computation " << cellList->computations << " times" <<endl;
    t1 = clock();
    cellList->computeCellList(Configuration->returnPositions());
    t2 = clock();
    cout << endl << "cellList took " << (t2-t1)/(scalar)CLOCKS_PER_SEC << " and iterated through the computation " << cellList->computations << " times" <<endl;

    ArrayHandle<dVec> p(Configuration->returnPositions());
    dVec target = p.data[0];
    printdVec(target);
    int cell = cellList->positionToCellIndex(target);
    ArrayHandle<unsigned int> particlesPerCell(cellList->elementsPerCell);
    ArrayHandle<int> indices(cellList->particleIndices);
    int neighs = particlesPerCell.data[cell];
    for (int nn = 0; nn < particlesPerCell.data[cell]; ++nn)
        {
        cout << "cell entry " << nn+1 << " out of "<< neighs << ": " << indices.data[cellList->cellListIndexer(nn,cell)] << endl;
        printdVec(p.data[indices.data[cellList->cellListIndexer(nn,cell)]]);
        };

    /*
    int cell;
    dVec bDims;
    PBC->getBoxDims(bDims);
    cell = cellList->positionToCellIndex(bDims);
    vector<int> cellNeighs;
    cellList->getCellNeighbors(cell,2,cellNeighs);
    for (int cc = 0; cc < cellNeighs.size(); ++cc)
        {
        iVec cellSet = cellList->indexToiVec(cellNeighs[cc]);
        cout <<cc<< ":  " << cellNeighs[cc] << endl;
        printiVec(cellSet);
        };

    bDims = 0.5*bDims;
    cell = cellList->positionToCellIndex(bDims);
    cellList->getCellNeighbors(cell,1,cellNeighs);
    for (int cc = 0; cc < cellNeighs.size(); ++cc)
        {
        iVec cellSet = cellList->indexToiVec(cellNeighs[cc]);
        cout <<cc<< ":  " << cellNeighs[cc] << endl;
        printiVec(cellSet);
        };
        */
/*
    IndexDD indexer(floor(L));

    for (int i1 =0; i1 < floor(L); ++i1)
        for (int i2 =0; i2 < floor(L); ++i2)
            for (int i3 =0; i3 < floor(L); ++i3)
                {
                iVec tester;
                tester.x[0] = i1;
                tester.x[1] = i2;
                tester.x[2] = i3;
                int result = indexer(tester) ;
                cout << tester.x[0] <<", "<< tester.x[1]<<", " << tester.x[2] <<" =   " << result << endl;
                cout << indexer.inverseIndex(result).x[0] << ", " << indexer.inverseIndex(result).x[1]  << ", " << indexer.inverseIndex(result).x[2] << endl;
                };
    cout << indexer.getNumElements() << "total elements" << endl;
*/
    /*
    iVec it(-1);it.x[0]-=1;
    iVec min(-1);
    iVec max(1);max.x[2]+=3;
    while(iVecIterate(it,min,max))
            printiVec(it);
    */

//
//The end of the tclap try
//
    } catch (ArgException &e)  // catch any exceptions
    { cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }
    return 0;
};
