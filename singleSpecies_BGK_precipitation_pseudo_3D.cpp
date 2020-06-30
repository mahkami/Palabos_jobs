/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2012 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* Main author: Babak, Chris and Andrea
 */

#include <cstdlib>
#include <iostream>

#include "palabos2D.h"
#include "palabos2D.hh"

using namespace plb;
using namespace std;

typedef double T;

#define NSDESCRIPTOR descriptors::ForcedD2Q9Descriptor
#define ADESCRIPTOR  descriptors::AdvectionDiffusionWithSourceD2Q9Descriptor

#define ADYNAMICS   AdvectionDiffusionWithSourceBGKdynamics
#define NSDYNAMICS GuoExternalForceBGKdynamics

#include "./porousMedia.h"
#include "./porousMedia.hh"
#include "./singleSpeciesHeterogeneousChemicalReactionsProcessingFunctional2D.h"
#include "./singleSpeciesHeterogeneousChemicalReactionsProcessingFunctional2D.hh"


std::string outDir, geometry;
T deltaRho;
T thickNess;
// T lx=735;
// T ly=555;
T lx=1680;
T ly=2180;

// we gp tp lattice units here
plint resolution ;
bool startsFromCheckingPoint;

struct Param {
    // dimensionless number ??
    
    T Pe; //Peclet number ?
    T Da; // Damköhler??
    plint poreSize;
    
    
    T trelaxNS; // relaxation time fluid 
    T trelaxSpecies1; // relaxation time primarie species
    T nu; //fluid viscosity
    
    plint maxIter;
    plint statIter;
    plint vtkIter;
    plint saveIter;
    
    plint parameterPorosity;
    plint nx,ny;
    
    // fluid stuff
    //T force ; //driving force for the fluid (based on reynolds definition???)
    T iniRho ;
    
    // BABAK CONSTANTS
    plint numberOfPrimariesSpecies; 
 
    T alpha ;
    T imposedConcentration1, initialConcentration1;
    T inletSpecies1; 
    T M;
    T Ke5;
    plint timeSplittingFluidReactions;
    
    Array<T,2>vectorForce;
    plint flowDirection;
    
    bool nonActiveNodesBounceBack;
    
} param;


void setParam(void)
{
        
    
    param.maxIter   = 10000001;
    param.statIter  =    1000;
    param.vtkIter   =     10000;
    param.saveIter  =    10000;
    
    
    param.trelaxNS = 1.;
    param.nu    = (param.trelaxNS-0.5)/NSDESCRIPTOR<T>::invCs2;

    
    param.vectorForce =  Array<T,2>(0.,0);
    param.flowDirection = 1;
    param.iniRho = 1.;
    
    param.numberOfPrimariesSpecies = 1;
        
    param.M = 1.;
    param.poreSize=45; //Width of main fracture in porous medium
       
    param.imposedConcentration1 = 0.1;
    param.initialConcentration1 = 0.1;

    param.inletSpecies1         = 1;
    
    param.timeSplittingFluidReactions =  1;
    
    param.nonActiveNodesBounceBack = true; // if you choose true you will get BB nodes for the advection diffusion
}


// /// A functional, used to initialize a pressure boundary to constant density
 template<typename T>
class ConstantDensity {
public:
    ConstantDensity(T density_)
        : density(density_)
    { }
    T operator()(plint iX, plint iY) const {
        return density;
    }
private:
    T density;
};

T computePermeability( MultiBlockLattice2D<T,NSDESCRIPTOR>& nsLattice, T nu , T deltaRho, plint flowDir, Box2D domain )
{
        //pcout << "Computing the permeability." << endl;

        // Compute only the y-direction of the velocity (direction of the flow).
        plint ny = nsLattice.getNy();

        T meanU = computeAverage(*computeVelocityComponent (nsLattice, domain, flowDir) );

        T k = (param.nu*meanU) / ((deltaRho/NSDESCRIPTOR<T>::invCs2)/(T)(ny-1)); 


        return k;
}

void problemSetupFluid( MultiBlockLattice2D<T,NSDESCRIPTOR>& nsLattice,
                      OnLatticeBoundaryCondition2D<T,NSDESCRIPTOR>* boundaryCondition, T deltaRho)
{
    
const plint nx = param.nx;
        const plint ny = param.ny;

        pcout << "Definition of inlet/outlet." << endl;
        Box2D top (0,nx-1, ny-1,ny-1);
        Box2D down(0,nx-1, 0,0);
        
	boundaryCondition->addPressureBoundary1P(top, nsLattice);
        boundaryCondition->addPressureBoundary1N(down, nsLattice);
	setBoundaryDensity (nsLattice, top,    ConstantDensity<T>(T(1)-deltaRho/2.) );
	setBoundaryDensity (nsLattice, down,   ConstantDensity<T>(T(1)+deltaRho/2.) );

        initializeAtEquilibrium( nsLattice, nsLattice.getBoundingBox(), 1., Array<T,2> (0.,0.)  );
    // the force has to be initialized
//    setExternalVector(nsLattice, nsLattice.getBoundingBox(),
//                      NSDESCRIPTOR<T>::ExternalField::forceBeginsAt, T deltaRho);
    
}

void setSpeciesAD(T darcyVel,  
                        MultiBlockLattice2D<T,ADESCRIPTOR>& mySpecies1,
                        MultiScalarField2D<T> porousMedium,
                        OnLatticeAdvectionDiffusionBoundaryCondition2D<T,ADESCRIPTOR>& species1BoundaryCondition)
{

    T D = std::fabs(darcyVel)*param.poreSize/param.Pe;   // poresize (45 in lattice units)
    
    param.alpha = D*param.Da*param.initialConcentration1/param.poreSize/param.poreSize;
    
    param.trelaxSpecies1 = D*ADESCRIPTOR<T>::invCs2 + 0.5 ;
    T omegaAD = 1./param.trelaxSpecies1;
    
    pcout << "New param.trelaxSpecies1= " << param.trelaxSpecies1 << std::endl;
    
    defineDynamics(mySpecies1, mySpecies1.getBoundingBox(), new ADYNAMICS<T, ADESCRIPTOR>(1./param.trelaxSpecies1) );
    
    initializeAtEquilibrium(mySpecies1, mySpecies1.getBoundingBox(), param.initialConcentration1, Array<T,2>(0.,0.) );
 
    
    Box2D left      (0, 0,  0,    param.ny-1);
    Box2D right     (param.nx-1, param.nx-1, 0,  param.ny-1);
    Box2D top       (0, param.nx-1,  param.ny-1,    param.ny-1);
    Box2D down      (0, param.nx-1,  0,   0);
    
    species1BoundaryCondition.addTemperatureBoundary1P(top,  mySpecies1);
    species1BoundaryCondition.addTemperatureBoundary1N(down, mySpecies1);
    
    setBoundaryDensity(mySpecies1, top,     param.initialConcentration1);
    setBoundaryDensity(mySpecies1, down,    param.inletSpecies1 );
    
        
//     if (param.nonActiveNodesBounceBack)
//         defineDynamics(mySpecies1, porousMedium, new BounceBack<T, ADESCRIPTOR>(param.initialConcentration1), 2 );
//     else
//         defineDynamics(mySpecies1, porousMedium, new NoDynamics<T, ADESCRIPTOR>(param.initialConcentration1), 2 );
    
    
    // initialize the dynamic of the adLattices based on flagPorousMedium
     initializeDynamicSinglepecies(mySpecies1, porousMedium, param.nonActiveNodesBounceBack, param.imposedConcentration1, omegaAD);
    

}


void writeVTKlattices(MultiBlockLattice2D<T,NSDESCRIPTOR>& nsLattice,
                        MultiBlockLattice2D<T,ADESCRIPTOR>& species1,
                        MultiScalarField2D<T> saturationIndex,
                        MultiScalarField2D<T> solidFraction,
                        MultiScalarField2D<T> changeInSolidFraction,
                        MultiScalarField2D<T> flagPorMedium,
                        MultiScalarField2D<T> porousMedium,
                        plint iter)
{
//     T dx = parameters.getDeltaX();
//     T dt = parameters.getDeltaT();
    
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), 1.);
    // Temperature is the order-0 moment of the advection-diffusion model. It can 
    //    therefore be computed with the function "computeDensity".
    vtkOut.writeData<float>(*computeDensity(species1), "concPrimarySpecies1", (T)1);
    vtkOut.writeData<2,float>(*computeVelocity(nsLattice), "velocity", 1.);
    vtkOut.writeData<float>(*computeDensity(nsLattice), "density", 1.);
    vtkOut.writeData<float>(flagPorMedium, "flagPorMedium", (T)1);
    vtkOut.writeData<float>(solidFraction, "solidFraction", (T)1);    

    
//     setToConstant(population, population.getBoundingBox(), 0.);
//     computePopulation(species1, population, species1.getBoundingBox(), 3);
    
//     std::string outfile_1 = global::directories().getOutputDir()+createFileName("concentration_", iter, 6)+".dat";       
//     plb_ofstream fout_1(outfile_1.c_str());
//     fout_1 << *computeDensity(species1)  << std::endl;
//     
//     std::string outfile_2 = global::directories().getOutputDir()+createFileName("solidFraction", iter, 6)+".dat";       
//     plb_ofstream fout_2(outfile_2.c_str());
//     fout_2 << *extractSubDomain(solidFraction, solidFraction.getBoundingBox())  << std::endl;
//         
//     std::string outfile_3 = global::directories().getOutputDir()+createFileName("flagPorMedium", iter, 6)+".dat";       
//     plb_ofstream fout_3(outfile_3.c_str());
//     fout_3 << *extractSubDomain(flagPorMedium, flagPorMedium.getBoundingBox())  << std::endl;
//         
//     std::string outfile4 = global::directories().getOutputDir()+createFileName("changeSolidFraction", iter, 6)+".dat";       
//     plb_ofstream fout4(outfile4.c_str());
//     fout4 << *extractSubDomain(changeInSolidFraction, changeInSolidFraction.getBoundingBox()) << std::endl;
//     
//     std::string outfile5 = global::directories().getOutputDir()+createFileName("normVel1_", iter, 6)+".dat";       
//     plb_ofstream fout5(outfile5.c_str());
//     fout5 << *computeVelocityNorm(nsLattice)  << std::endl;
//         
//     std::string outfile6 = global::directories().getOutputDir()+createFileName("porousMedium", iter, 6)+".dat";       
//     plb_ofstream fout6(outfile6.c_str());
//     fout6 << *extractSubDomain(porousMedium, porousMedium.getBoundingBox()) << std::endl;
// 	
// 	std::string outfile7 = global::directories().getOutputDir()+createFileName("density", iter, 6)+".dat";       
//     plb_ofstream fout7(outfile7.c_str());
//     fout7 << *computeDensity(nsLattice)  << std::endl;
    
}

void writeGif(MultiBlockLattice2D<T,NSDESCRIPTOR>& nsLattice,
              MultiBlockLattice2D<T,ADESCRIPTOR>& species1,int iT)

{
    const plint imSize = 600;
    const plint nx = nsLattice.getNx();
    const plint ny = nsLattice.getNy();
    Box2D slice(0, nx-1, 0, ny-1);
    ImageWriter<T> imageWriter("leeloo.map");
    imageWriter.writeScaledGif(createFileName("u", iT, 6),
                               *computeVelocityNorm(nsLattice, slice),
                               imSize, imSize);
    // Temperature is the order-0 moment of the advection-diffusion model. It can 
    //    therefore be computed with the function "computeDensity".
    imageWriter.writeScaledGif(createFileName("conc1ps", iT, 6),
                               *computeDensity(species1, slice),
                               imSize, imSize);
    
                         
                               
                               
                               
}

void writeGifsPM(MultiScalarField2D<T>& porousMedium,int iT)
{
    ImageWriter<T> imageWriter("leeloo.map");
    imageWriter.writeScaledGif(createFileName("porousMedium_", iT, 6), porousMedium);
}


int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    
    global::timer("simTime").start();
    srand (10);
    
    try {
        global::argv(1).read(outDir);
        global::argv(2).read(param.Pe);
        global::argv(3).read(param.Da);
        global::argv(4).read(deltaRho);
        global::argv(5).read(thickNess);
        global::argv(6).read(startsFromCheckingPoint);
        global::argv(7).read(geometry);
    }
    catch(PlbIOException& exception) {
        pcout << exception.what() << endl;
        pcout << "The structure of the input parameters should be : "
              << (string)global::argv(0) << " tmp Peclet (e.g. 10) Da (e.g. 0.5) deltaRho (e.g. 1e-3) thickness (e.g. 4) checkPoint 1 (1 if starting CheckPoint) geometry.txt" << endl;;
        // Exit the program, because wrong input data is a fatal error.
        exit(1);
    }
    
    
    param.nx = lx;
    param.ny = ly;
    
    
    std::string outputDirectory = outDir+"/";
    global::directories().setOutputDir(outputDirectory.c_str());

    std::string outfile_permeability = global::directories().getOutputDir()+"permeabilityEvolution.txt";       
    plb_ofstream fout_k(outfile_permeability.c_str());

    setParam();
           
    // Fluid Lattice 
    MultiBlockLattice2D<T, NSDESCRIPTOR> nsLattice (
            param.nx,param.ny,new NSDYNAMICS<T, NSDESCRIPTOR>(1./param.trelaxNS) );
    // Use periodic boundary conditions.
    nsLattice.periodicity().toggleAll(false);
            
    
    // LATTICES FOR THE PRIMARIES SPECIES ADVECTED BY THE FLUID FLOW
//     MultiBlockLattice2D<T, ADESCRIPTOR> species1 ( 
//             param.nx,param.ny, new NoDynamics<T, ADESCRIPTOR>(1./param.trelaxSpecies1) );
            
    MultiBlockLattice2D<T, ADESCRIPTOR> species1 ( 
            param.nx,param.ny, new ADYNAMICS<T, ADESCRIPTOR>(1./param.trelaxSpecies1) );
    // Use periodic boundary conditions.
    species1.periodicity().toggleAll(false);
                
    // We are using regularized boundary conditions 
    OnLatticeAdvectionDiffusionBoundaryCondition2D<T,ADESCRIPTOR>*
        species1BoundaryCondition = createLocalRegularizedAdvectionDiffusionBoundaryCondition2D<T,ADESCRIPTOR>();
    
    // initialization !!!!!!!!!!!!!!!!!!!!!!!!!!! for the fluid lattice.. 
    problemSetupFluid(nsLattice, createLocalBoundaryCondition2D<T,NSDESCRIPTOR>(), deltaRho);        
    // the species will be initialized once the steady state for the fluid is reached.        
    
    // MY REAL CHEMICAL SPECIES 
    MultiScalarField2D<T> saturationIndex(param.nx,param.ny); 
    saturationIndex.periodicity().toggleAll(true);
    
    setToConstant(saturationIndex, saturationIndex.getBoundingBox(), T());

    // Let's build the porous Medium (texture ScalarField),
    // 1st We bulid the porous matrix (texture field) by reading the sinthetic porous medium file
    // 2nd we set-up the dynamic for the fluid lattice (constructPorousMedia method) 
    // 3rd we initialize the flagMatrixes (active reaction node (BGKAdvDiffWithSource) (initialize flagStatus method) TO DO ! 
    // or solid (BB or Dirichelet Phase Field))
    // 4th properly i initialize the flagMatrixes (active reaction node (BGKAdvDiffWithSource) (initialize flagStatus method) TO DO !
    // 5th also the solidFraction can be initialized here
    
    MultiScalarField2D<T> texture(param.nx,param.ny);
    texture.periodicity().toggleAll(true);
    setToConstant(texture, texture.getBoundingBox(), 0.);   
    
    plb_ifstream textureFile(geometry.c_str());

    if (!textureFile.is_open()){
        pcout << "Error: could not open texture file " << endl;
    return -1;	
    }
    pcout << "Reading texture file" << endl;
    asciiStreamToUnSerializer<T>( &textureFile.getOriginalStream(),
    texture.getBlockUnSerializer(
    texture.getBoundingBox(), 
    global::IOpolicy().getIndexOrderingForStreams()
    ) );

    
    pcout << "Time spent for reading texture file: "
        << global::timer("texture").stop() << endl;    
        
    // Change the Dynamic of the fluid lattice from GuoExternalForceBGKdynamics to BounceBack where needed (on xtals nodes)
    constructPorousMedia(nsLattice, texture, 0);
    
    MultiScalarField2D<T> porousMedium(param.nx,param.ny);
    porousMedium.periodicity().toggleAll(true);
    
    perPorousMediaVisulatization(texture, porousMedium, 0);
    writeGifsPM(porousMedium,0);
   
    
    
    //////////////////  AND THEN I INITIALIZE THE SOLID FRACTION
    
  
    MultiScalarField2D<T> solidFraction(param.nx,param.ny); 
    solidFraction.periodicity().toggleAll(true);
    
    MultiScalarField2D<T> oldSolidFraction(param.nx,param.ny); 
    oldSolidFraction.periodicity().toggleAll(true);
    
    MultiScalarField2D<T> changeInSolidFraction(param.nx,param.ny); 
    changeInSolidFraction.periodicity().toggleAll(true);
    
    initializeSolidFraction(nsLattice,solidFraction);
    initializeSolidFraction(nsLattice,oldSolidFraction);
    
    setToConstant(changeInSolidFraction,changeInSolidFraction.getBoundingBox(),T());
    
    //////////////////  AND NOW THAT THE SOLID FRACTIOn IS WELL INITIALIZED, I INITIALIZE THE FLAGSTATUS 
    
    MultiScalarField2D<T> flagPorousMedium(param.nx,param.ny); // current time step nodes states 
    flagPorousMedium.periodicity().toggleAll(true);
    
    MultiScalarField2D<T> flagPorousMediumOld(param.nx,param.ny); // previous time step nodes states
    flagPorousMediumOld.periodicity().toggleAll(true);
    
    initializeFlagStatusTest(species1,solidFraction,flagPorousMedium);
    initializeFlagStatusTest(species1,solidFraction,flagPorousMediumOld);
    
    
    // Turn off internal statistics in order to improve parallel efficiency
    //   (it is not used anyway).
    nsLattice.toggleInternalStatistics(false);
    species1.toggleInternalStatistics(false);  
                        
    // for the chemical Species, you give an initial guess for primaries species and some other chemical species1
    // then you call the Speciation Processor. Once the Speciation is done you have all the chemical species that are well initialized
                        
                        
    Box2D everythingButInletAndOutlet (1,param.nx-2, 1,param.ny-2);
                        
    std::vector<MultiBlock2D* > blockPerDirichelet;
    blockPerDirichelet.push_back(&species1);
    blockPerDirichelet.push_back(&flagPorousMedium);
        
    
        std::vector<MultiBlock2D* > blockCalculateSourceSink;
    blockCalculateSourceSink.push_back(&species1);
    blockCalculateSourceSink.push_back(&solidFraction);
    blockCalculateSourceSink.push_back(&saturationIndex);
    blockCalculateSourceSink.push_back(&changeInSolidFraction);
    blockCalculateSourceSink.push_back(&flagPorousMedium);
    
        std::vector<MultiBlock2D* > blockUpdateFlag;
    blockUpdateFlag.push_back(&species1);
    blockUpdateFlag.push_back(&solidFraction);
    blockUpdateFlag.push_back(&changeInSolidFraction);
    blockUpdateFlag.push_back(&flagPorousMedium);
    blockUpdateFlag.push_back(&flagPorousMediumOld);
    blockUpdateFlag.push_back(&nsLattice);
    
    std::vector<MultiBlock2D* > blockUpdateSpecies;
    blockUpdateSpecies.push_back(&species1);
    blockUpdateSpecies.push_back(&changeInSolidFraction);
    blockUpdateSpecies.push_back(&flagPorousMedium);
    blockUpdateSpecies.push_back(&flagPorousMediumOld);
    
    
    // fluid related shit
    std::vector<MultiBlock2D* > blockChangeDynamicFluid;
    blockChangeDynamicFluid.push_back(&nsLattice);
    blockChangeDynamicFluid.push_back(&species1);
    blockChangeDynamicFluid.push_back(&solidFraction);
    blockChangeDynamicFluid.push_back(&flagPorousMedium);
    blockChangeDynamicFluid.push_back(&flagPorousMediumOld);
    
    
    std::vector<MultiBlock2D* > blockUpdVelChemSpecies;
    blockUpdVelChemSpecies.push_back(&nsLattice);
    blockUpdVelChemSpecies.push_back(&species1);


    util::ValueTracer<T> steadyStateReached((T)1,(T)10, 1e-4);
    bool convergedOnce = false;

    pcout << "Porosity: " << computePorosity(solidFraction, 0.5, solidFraction.getBoundingBox()) << endl;    

    
    plint iniTime = 0;
    
    if (startsFromCheckingPoint)
    {
        loadBinaryBlock(nsLattice, "checkpoint.dat");
        pcout << "We start from checkingPoint for velocity "  << iniTime << " !!! We are ready to rockify." << endl;
        convergedOnce = true;
        T darcyVelocityConverged = 
                        computeDarcyVelocity(nsLattice, texture, Box2D(0,param.nx-1,param.ny/2,param.ny/2),
                                            param.parameterPorosity, param.flowDirection);
           
        pcout << "darcyVelocityConverged: " << darcyVelocityConverged <<  std::endl;
       
        T porosity = computePorosity(solidFraction, 0.5, solidFraction.getBoundingBox());
            pcout << "New Porosity: " << porosity << std::endl;
        
        T permeability = computePermeability(nsLattice, param.nu, deltaRho, param.flowDirection, nsLattice.getBoundingBox());
        fout_k << iniTime  << " " << porosity <<  " " << permeability << std::endl;
        setSpeciesAD(darcyVelocityConverged, species1, flagPorousMedium,*species1BoundaryCondition);
        species1.initialize();   
       
        writeVTKlattices(nsLattice, species1, saturationIndex, solidFraction,
                            changeInSolidFraction, flagPorousMedium, porousMedium, iniTime);
        pcout << "Writing gif..." << endl;
        writeGif(nsLattice,species1,iniTime);
   
    }
    
    
    
    
    // Main loop over time iterations.
    for (plint iT = 0; iT < param.maxIter; ++iT) 
    {
        
//        
        
         if (iT % param.saveIter == 0  && convergedOnce == false  ) {
            // dissolution and reaction processes are not active at this point
                      
            T darcyVelocity = computeDarcyVelocity(nsLattice, texture, Box2D(0,param.nx-1,param.ny/2,param.ny/2), param.parameterPorosity, param.flowDirection );
            pcout << "darcyVel " << darcyVelocity << std::endl;
            //writeGif(nsLattice,species1,iT);
            steadyStateReached.takeValue(darcyVelocity,true);
        }
        
        
        if(steadyStateReached.hasConverged() && convergedOnce == false  )
        {
            
            iniTime = iT;
            convergedOnce = true;
            saveBinaryBlock(nsLattice, "checkpoint.dat"); 
            
            T darcyVelocityConverged =  computeDarcyVelocity(nsLattice, texture, Box2D(0,param.nx-1,param.ny/2,param.ny/2), param.parameterPorosity, param.flowDirection);
            
            pcout << "darcyVelocityConverged: " << darcyVelocityConverged <<  std::endl; 
            
            
            T porosity = computePorosity(solidFraction, 0.5, solidFraction.getBoundingBox());
            pcout << "New Porosity: " << porosity << std::endl;
            
            T permeability = (computePermeability(nsLattice, param.nu, deltaRho, param.flowDirection, nsLattice.getBoundingBox()));
            fout_k << iT << " " << porosity << " " << permeability << std::endl;
			
            pcout << "We have reached the steady state (Darcy) at iteration "  << iT << " !!! We are ready to rockify." << endl;
            
            // why are we using darcyVelocityConverged? instead of plain darcy velocity? values different significantly.
            setSpeciesAD(darcyVelocityConverged, species1, flagPorousMedium,*species1BoundaryCondition);
            species1.initialize();    
            
            writeVTKlattices(nsLattice, species1, 
                             saturationIndex, solidFraction, changeInSolidFraction, flagPorousMedium, porousMedium, iT-iniTime-1);
            pcout << "Writing gif..." << endl;
            // writeGif(nsLattice,species1,iT-iniTime-1);
            
        }


        if (iT % param.saveIter == 0 && convergedOnce )
        {
            pcout << "At iter " << iT  << std::endl;
            pcout << "Writing VTK..." << endl;
            porousMediaEvolution(nsLattice,porousMedium);
            writeVTKlattices(nsLattice, species1, 
                             saturationIndex, solidFraction, changeInSolidFraction, flagPorousMedium, porousMedium, iT-iniTime);
            pcout << "Writing gif..." << endl;
            // writeGif(nsLattice,species1,iT-iniTime);
            T porosity = computePorosity(solidFraction, 0.5, solidFraction.getBoundingBox());
            pcout << "New Porosity: " << porosity << std::endl;
            
            T permeability = computePermeability (nsLattice, param.nu, deltaRho, param.flowDirection, nsLattice.getBoundingBox());
            fout_k << iT  << " " << porosity <<  " " << permeability << std::endl;
        }
        
        if (convergedOnce)
        {
            
            species1.collideAndStream();
            
            if (param.nonActiveNodesBounceBack == false)
            {
                applyProcessingFunctional (
                            new ChrisDiricheletProcessor2D<T,ADESCRIPTOR> (param.imposedConcentration1, 1./param.trelaxSpecies1),
                            everythingButInletAndOutlet,
                            blockPerDirichelet);    
            }
            
            applyProcessingFunctional (
                new CalculateSourceSinkAndUpdateSolidFractionProcessor2D<T,ADESCRIPTOR> (param.Ke5, param.alpha, param.M, param.imposedConcentration1),
                everythingButInletAndOutlet,
                blockCalculateSourceSink);    
                
            applyProcessingFunctional (
                new UpdateFlagProcess2D<T,ADESCRIPTOR,NSDESCRIPTOR> (),
                everythingButInletAndOutlet,
                blockUpdateFlag);
    // 
            applyProcessingFunctional (
                new CalculateSourceSinkPerPrimariesSpeciesProcessor2D<T,ADESCRIPTOR> (param.M, 1./param.trelaxSpecies1, 
                                                                                    param.imposedConcentration1, param.nonActiveNodesBounceBack),
                everythingButInletAndOutlet,
                blockUpdateSpecies);

            applyProcessingFunctional (
                    new FlatAdiabaticBoundaryFunctional2D<T,ADESCRIPTOR,1,1> (),
                        Box2D(0,param.nx-1,param.ny-1,param.ny-1), species1);  


            
        }            
        nsLattice.collideAndStream();
        
       
        
                
        applyProcessingFunctional (
            new ChangeDynamicDissolutionPrecipitationProcessor2D<T,NSDESCRIPTOR,ADESCRIPTOR> (1./param.trelaxNS, param.vectorForce, param.iniRho),
            everythingButInletAndOutlet,
            blockChangeDynamicFluid);
        
        
         applyProcessingFunctional(   
            new DragTerm2D<T,NSDESCRIPTOR> (thickNess, param.nu),
                everythingButInletAndOutlet, nsLattice);
        
        

//     // WE SEND THE NEW VELOCITIES TO THE PRIMARIES SPECIES, READY TO BE ADVECTED AGAIn        
        applyProcessingFunctional (
            new PassiveScalarProcessor2D<T,NSDESCRIPTOR,ADESCRIPTOR> (),
            nsLattice.getBoundingBox(),
            blockUpdVelChemSpecies);
    

    }
    
    writeGif(nsLattice,species1,param.maxIter);
    
//     T tEnd = global::timer("simTime").stop();
    

}
