// Authors: 
// Christian Huber, chris_huber@berkeley.edu
// Jonas Latt, jonas@lbmethod.org
// Andrea Parmigiani, Andrea.Parmigiani@unige.ch
//
// Copyright (c) 2009 by Ecole Polytechnique Federale de Lausanne (EPFL)
// 1015 Lausanne
// Switzerland
//
// All rights reserved.

#ifndef POROUS_MEDIA_HH
#define POROUS_MEDIA_HH

#include "porousMedia.h"
// #include "./chemicalReactionsProcessors/chrisDiricheletBoundaryDynamics.h"
// #include "./chemicalReactionsProcessors/chrisDiricheletBoundaryDynamics.hh"

template<typename T, template<typename U> class Descriptor>
ConstructPorousMediaFunctional2D<T,Descriptor>::ConstructPorousMediaFunctional2D (plint parameterPorosity_ )
    : parameterPorosity(parameterPorosity_)
{ }

template<typename T, template<typename U> class Descriptor>
void ConstructPorousMediaFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& geometry )
{
    Dot2D relativeOffset = computeRelativeDisplacement(lattice,geometry);
    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint oX = relativeOffset.x + iX;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint oY = relativeOffset.y + iY;
                plint materialIndex = util::roundToInt(geometry.get(oX,oY));
                //if (materialIndex >= parameterPorosity) {  // Solid cell    for geometryBabak
                if (materialIndex <= parameterPorosity) {  // Solid cell   // for newTexture    
                    //lattice.attributeDynamics(iX,iY, new BounceBack<T,Descriptor>(T(1)));
                    defineDynamics(lattice, iX,iY,new BounceBack<T,Descriptor>()); 
                }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
ConstructPorousMediaFunctional2D<T,Descriptor>* ConstructPorousMediaFunctional2D<T,Descriptor>::clone() const {
    return new ConstructPorousMediaFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ConstructPorousMediaFunctional2D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables;
    modified[1] = modif::nothing;
}


template<typename T, template<typename U> class Descriptor>
void constructPorousMedia(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& texture, plint parameterPorosity)
{
    applyProcessingFunctional (
        new ConstructPorousMediaFunctional2D<T,Descriptor>(parameterPorosity),
        lattice.getBoundingBox(),
        lattice, texture );
}

/////////////////// INITIALIZE FLAG STATUS


template<typename T, template<typename U> class Descriptor>
void InitializeFlagStatusFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& flagStatus )
{
   
    Dot2D relativeOffsetPM = computeRelativeDisplacement(lattice,flagStatus);
  
    BounceBack<T,Descriptor> BBdynamics;
    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) 
    {
        plint pmX = relativeOffsetPM.x + iX;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) 
        {
            plint pmY = relativeOffsetPM.y + iY;
            
            Cell<T,Descriptor>& cellNS  = lattice.get(iX,iY);
            
            if( cellNS.getDynamics().getId() != BBdynamics.getId()  )
            {
                flagStatus.get(pmX,pmY) = T(0); // 0 stands for liquid Node;
                
            }
            else if ( cellNS.getDynamics().getId() == BBdynamics.getId()  )
            {   
                plint activeLatticeNode = 0;
                for (plint iPop = 1; iPop< Descriptor<T>::q; iPop++)
                {
                    Cell<T,Descriptor>& cellNeigh = lattice.get(iX+Descriptor<T>::c[iPop][0],iY + Descriptor<T>::c[iPop][1]);
                    if (cellNeigh.getDynamics().getId() != BBdynamics.getId() )
                    {
                        activeLatticeNode += 1; // 0 stands for Solid but active node (where BGK with sink source term is still active);
                    }
                }
                
                if (activeLatticeNode >= 1 )
                    flagStatus.get(pmX,pmY) = T(1); // 0 stands for Solid but active node (where BGK with sink source term is still active);
                else
                    flagStatus.get(pmX,pmY) = T(2);
                    
            
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
InitializeFlagStatusFunctional2D<T,Descriptor>* InitializeFlagStatusFunctional2D<T,Descriptor>::clone() const {
    return new InitializeFlagStatusFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void InitializeFlagStatusFunctional2D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}



template<typename T, template<typename U> class Descriptor>
void initializeFlagStatus(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& flagStatus)
{
    applyProcessingFunctional (
        new InitializeFlagStatusFunctional2D<T,Descriptor>(),
        lattice.getBoundingBox(),
        lattice, flagStatus );
}


///////  Cleaning porous media geometry



// template<typename T, template<typename U> class Descriptor>
// ConstructPorousMediaFunctional2D<T,Descriptor>::ConstructPorousMediaFunctional2D (plint parameterPorosity_ )
//     : parameterPorosity(parameterPorosity_)
// { }

template<typename T, template<typename U> class Descriptor>
CleaningFlagStatusFunctional2D<T,Descriptor>::CleaningFlagStatusFunctional2D(T iniRho_, Array<T,2> u_, T omega_ )
    : iniRho(iniRho_), u(u_), omega(omega_)
{ }


template<typename T, template<typename U> class Descriptor>
void CleaningFlagStatusFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& flagStatus )
{
   
    Dot2D relativeOffsetPM = computeRelativeDisplacement(lattice,flagStatus);
  
    BounceBack<T,Descriptor> BBdynamics;
    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) 
    {
        plint pmX = relativeOffsetPM.x + iX;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) 
        {
            plint pmY = relativeOffsetPM.y + iY;
            
            Cell<T,Descriptor>& cellNS  = lattice.get(iX,iY);
            
            if (flagStatus.get(pmX,pmY) == T(1))
            {
                bool myStatusOk = false; // I am a true active node
                for (plint iPop = 1; iPop< Descriptor<T>::q; iPop++)
                {
                    plint next_pmX = pmX + Descriptor<T>::c[iPop][0]; 
                    plint next_pmY = pmY + Descriptor<T>::c[iPop][1];
                    
                    //Cell<T,Descriptor>& cellNeigh = lattice.get(iX+Descriptor<T>::c[iPop][0],iY + Descriptor<T>::c[iPop][1]);
                    
                    if (flagStatus.get(next_pmX,next_pmY) == T(2) && myStatusOk == false )
                    {
                        myStatusOk = true; // 0 stands for Solid but active node (where BGK with sink source term is still active);
                    }
                    
                }
                
                if (myStatusOk == false)
                {
                    flagStatus.get(pmX,pmY) = T(0);
                    if (cellNS.getDynamics().getId() == BBdynamics.getId() )
                    {
                        defineDynamics(lattice,iX,iY, new RegularizedBGKdynamics<T,Descriptor>(omega));  
                        
                    
                        iniCellAtEquilibrium(cellNS, iniRho, u);

                        
                    }
                    
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
CleaningFlagStatusFunctional2D<T,Descriptor>* CleaningFlagStatusFunctional2D<T,Descriptor>::clone() const {
    return new CleaningFlagStatusFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void CleaningFlagStatusFunctional2D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables;
    modified[1] = modif::staticVariables;
}





template<typename T, template<typename U> class Descriptor>
void cleanFlagStatus(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& flagStatus, plint iniRho, Array<T,2> u, T omega)
{
    applyProcessingFunctional (
        new CleaningFlagStatusFunctional2D<T,Descriptor>(iniRho,u,omega),
        lattice.getBoundingBox(),
        lattice, flagStatus );
}






template<typename T>
ComputeSigmaValueFunctional2D<T>::ComputeSigmaValueFunctional2D(plint flagValue_)
    : numInterfaceCellId(this->getStatistics().subscribeIntSum()), 
    flagValue(flagValue_)
{ }

template<typename T>
void ComputeSigmaValueFunctional2D<T>::process (
        Box2D domain, ScalarField2D<T>& flagStatus )
{
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint materialIndex = util::roundToInt(flagStatus.get(iX,iY));
                if (materialIndex == flagValue) {  // Fluid Cell
                    statistics.gatherIntSum(numInterfaceCellId, 1);
                }
        }
    }
}

template<typename T>
ComputeSigmaValueFunctional2D<T>* ComputeSigmaValueFunctional2D<T>::clone() const
{
    return new ComputeSigmaValueFunctional2D<T>(*this);
}

template<typename T>
plint ComputeSigmaValueFunctional2D<T>::getNumInterfaceCells() const {
    return this->getStatistics().getIntSum(numInterfaceCellId);
}

template<typename T>
void ComputeSigmaValueFunctional2D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;
}





template<typename T>
T getSigmaValue(MultiScalarField2D<T>& flagStatus, Box2D domain, plint flagValue)
{
   
    ComputeSigmaValueFunctional2D<T> functional(flagValue);
    applyProcessingFunctional(functional, domain, flagStatus);
    
    T totalInterface = functional.getNumInterfaceCells();
    T totalCell = (domain.x1-domain.x0)*(domain.y1-domain.y0);
    
    return totalInterface / totalCell;
        
}



///////////////////////

template<typename T, template<typename U> class Descriptor>
void InitializeFlagStatusFunctionalTest2D<T,Descriptor>::
                processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks)
{
    
    
    BlockLattice2D<T,Descriptor>& species1 =
            *dynamic_cast<BlockLattice2D<T,Descriptor> *>(atomicBlocks[0]);
                   
    ScalarField2D<T>& solidFraction =
            *dynamic_cast<ScalarField2D<T> *>(atomicBlocks[1]);
           
    ScalarField2D<T>& flagStatus =
            *dynamic_cast<ScalarField2D<T> *>(atomicBlocks[2]);
            
            
    Dot2D relativeOffsetSF = computeRelativeDisplacement(species1,solidFraction);
    Dot2D relativeOffsetPM = computeRelativeDisplacement(species1,flagStatus);
  
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) 
    {
        plint sfX = relativeOffsetSF.x + iX;
        plint pmX = relativeOffsetPM.x + iX;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) 
        {
            plint sfY = relativeOffsetSF.y + iY;
            plint pmY = relativeOffsetPM.y + iY;
            
            
            
         
            T liquidFraction = T(1)-solidFraction.get(sfX,sfY);
            plint test = 0;
            plint tbb  = 0;
            
            for (plint iPop = 1; iPop< Descriptor<T>::q; iPop++)
            {

                T neighLiquidFraction = T(1)-solidFraction.get(sfX+Descriptor<T>::c[iPop][0],sfY+ Descriptor<T>::c[iPop][1]);
            
                if( liquidFraction < (T)1 && neighLiquidFraction == (T)1)
                    test=1; // Source Sink Active
                if( liquidFraction > (T)0 && neighLiquidFraction == (T)0 )
                    test=1; // Source Sink active
                if (liquidFraction < 0.5 && neighLiquidFraction < 0.5) //??
                    tbb += 1;
                
            }    
            
            if(tbb == Descriptor<T>::q - 1 ) // mines rest velocity 
                test=2;     
            
            flagStatus.get(pmX,pmY) = T(test); // 0 stands for Solid but active node (where BGK with sink source term is still active);
                    
            
        }
    }
}

template<typename T, template<typename U> class Descriptor>
InitializeFlagStatusFunctionalTest2D<T,Descriptor>* InitializeFlagStatusFunctionalTest2D<T,Descriptor>::clone() const {
    return new InitializeFlagStatusFunctionalTest2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void InitializeFlagStatusFunctionalTest2D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}


template<typename T, template<typename U> class Descriptor>
void initializeFlagStatusTest(MultiBlockLattice2D<T,Descriptor>& species1, MultiScalarField2D<T>& solidFraction, MultiScalarField2D<T>& flagStatus)
{
    std::vector<MultiBlock2D* > block;
    block.push_back(&species1);
    block.push_back(&solidFraction);
    block.push_back(&flagStatus);
    
    applyProcessingFunctional (
        new InitializeFlagStatusFunctionalTest2D<T,Descriptor>(),
        species1.getBoundingBox(),
        block );
}







// INITALIZE SOLID FRACTION


template<typename T, template<typename U> class Descriptor>
void InitializeSolifFractionFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& solidFraction )
{
   
    Dot2D relativeOffsetPM = computeRelativeDisplacement(lattice,solidFraction);
  
    BounceBack<T,Descriptor> BBdynamics;
    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) 
    {
        plint pmX = relativeOffsetPM.x + iX;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) 
        {
            plint pmY = relativeOffsetPM.y + iY;
            
            Cell<T,Descriptor>& cellNS  = lattice.get(iX,iY);
            
            if( cellNS.getDynamics().getId() == BBdynamics.getId()  )
                solidFraction.get(pmX,pmY) = T(1); // 
            else
                solidFraction.get(pmX,pmY) = T(0); // 
         
        }
    }
}

template<typename T, template<typename U> class Descriptor>
InitializeSolifFractionFunctional2D<T,Descriptor>* InitializeSolifFractionFunctional2D<T,Descriptor>::clone() const {
    return new InitializeSolifFractionFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void InitializeSolifFractionFunctional2D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}



template<typename T, template<typename U> class Descriptor>
void initializeSolidFraction(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& solidFraction)
{
    applyProcessingFunctional (
        new InitializeSolifFractionFunctional2D<T,Descriptor>(),
        lattice.getBoundingBox(),
        lattice, solidFraction );
}

//////////////////////  POROUS MEDIA VISUALIZATION

template<typename T>
VisualizePorousMedium2D<T>::VisualizePorousMedium2D (plint parameterPorosity_ )
    : parameterPorosity(parameterPorosity_)
{ }

template<typename T>
void VisualizePorousMedium2D<T>::process (
        Box2D domain, ScalarField2D<T>& texture, ScalarField2D<T>& porousMedium )
{
   
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) 
    {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) 
        {
            plint materialIndex = util::roundToInt(texture.get(iX,iY));
            //if (materialIndex >= parameterPorosity) //forGeometrybabak 
            if (materialIndex <= parameterPorosity) // for newTexture
            {
                porousMedium.get(iX,iY) = T(0); //solid
            }
            else
            {
                porousMedium.get(iX,iY) = T(1); //fluid
            }
        }
    }
}

template<typename T>
VisualizePorousMedium2D<T>* VisualizePorousMedium2D<T>::clone() const {
    return new VisualizePorousMedium2D<T>(*this);
}

template<typename T>
void VisualizePorousMedium2D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}



template<typename T>
void perPorousMediaVisulatization(MultiScalarField2D<T>& texture, MultiScalarField2D<T>& porousMedium, plint parameterPorosity)
{
    applyProcessingFunctional ( new VisualizePorousMedium2D<T>(parameterPorosity), texture.getBoundingBox(), texture, porousMedium );
}




template<typename T, template<typename U> class Descriptor>
InitializeDynamicsSpecies2D<T,Descriptor>::InitializeDynamicsSpecies2D (bool nonActiveNodesBounceBack_, T iniConcentration1_, 
                                                                        T iniConcentration2_, T iniConcentration3_)
    : nonActiveNodesBounceBack(nonActiveNodesBounceBack_), iniConcentration1(iniConcentration1_),
        iniConcentration2(iniConcentration2_), iniConcentration3(iniConcentration3_)
{ }



template<typename T, template<typename U> class Descriptor>
void InitializeDynamicsSpecies2D<T,Descriptor>::
        processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks)
{
 
 
    
    BlockLattice2D<T,Descriptor>& species1 =
            *dynamic_cast<BlockLattice2D<T,Descriptor> *>(atomicBlocks[0]);
            
    BlockLattice2D<T,Descriptor>& species2 =
            *dynamic_cast<BlockLattice2D<T,Descriptor> *>(atomicBlocks[1]);
            
    BlockLattice2D<T,Descriptor>& species3 =
            *dynamic_cast<BlockLattice2D<T,Descriptor> *>(atomicBlocks[2]);
            
    ScalarField2D<T>& porousMediumFlag =
            *dynamic_cast<ScalarField2D<T> *>(atomicBlocks[3]);
            
    
    
    Dot2D relativeOffsetPM = computeRelativeDisplacement(species1,porousMediumFlag);
  
 
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) 
    {
        plint pmX = relativeOffsetPM.x + iX;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) 
        {
            plint pmY = relativeOffsetPM.y + iY;
         
           
            if (util::roundToInt(porousMediumFlag.get(pmX,pmY)) == 2){
    
                if (nonActiveNodesBounceBack){
                    defineDynamics(species1, iX,iY, new BounceBack<T,Descriptor>(iniConcentration1));  
                    defineDynamics(species2, iX,iY, new BounceBack<T,Descriptor>(iniConcentration2));  
                    defineDynamics(species3, iX,iY, new BounceBack<T,Descriptor>(iniConcentration3)); 
//                     pcout << "do I enter here ?" << std::endl;
                }
                else{ // if the DiricheletPhaseFieldBoundary is chosen 
                    
                    defineDynamics(species1, iX,iY, new NoDynamics<T,Descriptor>(iniConcentration1));  
                    defineDynamics(species2, iX,iY, new NoDynamics<T,Descriptor>(iniConcentration2));  
                    defineDynamics(species3, iX,iY, new NoDynamics<T,Descriptor>(iniConcentration3)); 
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
InitializeDynamicsSpecies2D<T,Descriptor>*
    InitializeDynamicsSpecies2D<T,Descriptor>::clone() const
{
    return new InitializeDynamicsSpecies2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void InitializeDynamicsSpecies2D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    
    modified[0] = modif::dataStructure; // species1
    modified[1] = modif::dataStructure; // species2
    modified[2] = modif::dataStructure; // species3
    modified[3] = modif::nothing; // flagMedia
   
}




template<typename T, template<typename U> class Descriptor>
void initializeDynamicPrincipleSpecies(MultiBlockLattice2D<T,Descriptor>& species3, MultiBlockLattice2D<T,Descriptor>& species2,
                                       MultiBlockLattice2D<T,Descriptor>& species1, MultiScalarField2D<T>& porousMediumFlag,
                                       bool nonActiveNodesBounceBack, T iniConcentration1, T iniConcentration2, T iniConcentration3)
{

    std::vector<MultiBlock2D* > blockIniDynamic;
    blockIniDynamic.push_back(&species1);
    blockIniDynamic.push_back(&species2);
    blockIniDynamic.push_back(&species3);
    blockIniDynamic.push_back(&porousMediumFlag);
    
    
    
    applyProcessingFunctional (
            new InitializeDynamicsSpecies2D<T,Descriptor>(nonActiveNodesBounceBack, iniConcentration1, iniConcentration2, iniConcentration3),
                               species3.getBoundingBox(), blockIniDynamic );

}


//////////////////////////////////////////////



template<typename T, template<typename U> class Descriptor>
InitializeDynamicsSingleSpecies2D<T,Descriptor>::InitializeDynamicsSingleSpecies2D (bool nonActiveNodesBounceBack_, T imposedConcentration_, T omegaAD_)
    : nonActiveNodesBounceBack(nonActiveNodesBounceBack_), imposedConcentration(imposedConcentration_), omegaAD(omegaAD_)
{ }


template<typename T, template<typename U> class Descriptor>
void InitializeDynamicsSingleSpecies2D<T,Descriptor>::
        processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks)
{
 
 
    
    BlockLattice2D<T,Descriptor>& species1 =
            *dynamic_cast<BlockLattice2D<T,Descriptor> *>(atomicBlocks[0]);
             
            
    ScalarField2D<T>& porousMediumFlag =
            *dynamic_cast<ScalarField2D<T> *>(atomicBlocks[1]);
            
    
    Dot2D relativeOffsetPM = computeRelativeDisplacement(species1,porousMediumFlag);
  
                    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) 
    {
        plint pmX = relativeOffsetPM.x + iX;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) 
        {
            plint pmY = relativeOffsetPM.y + iY;
             
             if (porousMediumFlag.get(pmX,pmY) == 2){    
                if (nonActiveNodesBounceBack){
                    defineDynamics( species1, iX,iY, new BounceBack<T,Descriptor>(T(imposedConcentration)) );
                }
                else{ // if the DiricheletPhaseFieldBoundary is chosen 
                    
                    defineDynamics( species1, iX,iY, new NoDynamics<T,Descriptor>(T(imposedConcentration)) );
                
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
InitializeDynamicsSingleSpecies2D<T,Descriptor>*
    InitializeDynamicsSingleSpecies2D<T,Descriptor>::clone() const
{
    return new InitializeDynamicsSingleSpecies2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void InitializeDynamicsSingleSpecies2D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    
    modified[0] = modif::dynamicVariables; // species1
    modified[1] = modif::nothing; // flagMedia
   
}









template<typename T, template<typename U> class Descriptor>
void initializeDynamicSinglepecies(MultiBlockLattice2D<T,Descriptor>& species1, MultiScalarField2D<T>& porousMediumFlag,
                                       bool nonActiveNodesBounceBack, T imposedConcentration, T omegaAD)
{

    std::vector<MultiBlock2D* > blockIniDynamic;
    blockIniDynamic.push_back(&species1);
    blockIniDynamic.push_back(&porousMediumFlag);
    
    
    
    applyProcessingFunctional (
            new InitializeDynamicsSingleSpecies2D<T,Descriptor>(nonActiveNodesBounceBack,imposedConcentration, omegaAD),species1.getBoundingBox(), blockIniDynamic );

}



//////////////// Porous Media Evolution


template<typename T, template<typename U> class Descriptor>
void PorousMediaEvolution2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& porousMediumFlag )
{
 
            
    
    Dot2D relativeOffsetPM = computeRelativeDisplacement(lattice,porousMediumFlag);
  
  
    BounceBack<T,Descriptor> BBdynamics;
    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) 
    {
        plint pmX = relativeOffsetPM.x + iX;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) 
        {
            plint pmY = relativeOffsetPM.y + iY;
         
    
            Cell<T,Descriptor>&   cellNS = lattice.get(iX,iY);
            
            if (cellNS.getDynamics().getId() != BBdynamics.getId())
            {
                porousMediumFlag.get(pmX,pmY) = T(1);
            }
            else
            {
                porousMediumFlag.get(pmX,pmY) = T(0);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
PorousMediaEvolution2D<T,Descriptor>*
    PorousMediaEvolution2D<T,Descriptor>::clone() const
{
    return new PorousMediaEvolution2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void PorousMediaEvolution2D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    
    modified[0] = modif::nothing; // species1
    modified[1] = modif::staticVariables; // flagMedia
   
}


template<typename T, template<typename U> class Descriptor>
void porousMediaEvolution(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& porousMediumFlag)
{
     applyProcessingFunctional (
            new PorousMediaEvolution2D<T,Descriptor>(),lattice.getBoundingBox(),
        lattice, porousMediumFlag );

}


template<typename T>
ComputePorosityFunctional2D<T>::ComputePorosityFunctional2D(T parameterPorosity_)
    : numFluidCellId(this->getStatistics().subscribeIntSum()),
      parameterPorosity(parameterPorosity_)
{ }

template<typename T>
void ComputePorosityFunctional2D<T>::process (
        Box2D domain, ScalarField2D<T>& geometry )
{
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
//                 plint materialIndex = util::roundToInt(geometry.get(iX,iY));
                plint materialIndex = geometry.get(iX,iY);
                if (materialIndex <= parameterPorosity) {  // Fluid Cell
                //if (materialIndex < parameterPorosity) {  // Fluid Cell per babak geometry
                    statistics.gatherIntSum(numFluidCellId, 1);
                }
        }
    }
}

template<typename T>
ComputePorosityFunctional2D<T>* ComputePorosityFunctional2D<T>::clone() const
{
    return new ComputePorosityFunctional2D<T>(*this);
}

template<typename T>
plint ComputePorosityFunctional2D<T>::getNumFluidCells() const {
    return this->getStatistics().getIntSum(numFluidCellId);
}

template<typename T>
void ComputePorosityFunctional2D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;
}





///////////////////////

template<typename T, template<typename U> class Descriptor>
ComputePorosityFromFluidFunctional2D<T,Descriptor>::ComputePorosityFromFluidFunctional2D(bool flag_)
    : numFluidCellId(this->getStatistics().subscribeIntSum()),
      flag(flag_)
{ }



template<typename T, template<typename U> class Descriptor>
void ComputePorosityFromFluidFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    BlockStatistics& statistics = this->getStatistics();
    BounceBack<T,Descriptor> BBdynamics;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                Cell<T,Descriptor>&   cellNS = lattice.get(iX,iY);
                
                if (cellNS.getDynamics().getId() != BBdynamics.getId())
                {
                    statistics.gatherIntSum(numFluidCellId, 1);
                }
                
        }
    }
}

template<typename T, template<typename U> class Descriptor>
ComputePorosityFromFluidFunctional2D<T,Descriptor>* ComputePorosityFromFluidFunctional2D<T,Descriptor>::clone() const
{
    return new ComputePorosityFromFluidFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
plint ComputePorosityFromFluidFunctional2D<T,Descriptor>::getNumFluidCells() const {
    return this->getStatistics().getIntSum(numFluidCellId);
}

template<typename T, template<typename U> class Descriptor>
void ComputePorosityFromFluidFunctional2D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;
}

////////////////

template<typename T, template<typename U> class Descriptor>
T computePorosityFromFluidLattice(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain) {
    ComputePorosityFromFluidFunctional2D<T,Descriptor> functional(true);
    applyProcessingFunctional(functional, domain, lattice);
    return (T) functional.getNumFluidCells() / (T) domain.nCells();
}







template<typename T>
T computePorosity(MultiScalarField2D<T>& geometry, T parameterPorosity, Box2D domain) {
    ComputePorosityFunctional2D<T> functional(parameterPorosity);
    applyProcessingFunctional(functional, domain, geometry);
    return (T) functional.getNumFluidCells() / (T) domain.nCells();
}


template<typename T, template<typename U> class Descriptor>
DarcyVelocityFunctional2D<T,Descriptor>::DarcyVelocityFunctional2D(plint parameterPorosity_, plint flowDirection_)
    : sumDensityId   (this->getStatistics().subscribeSum()),
      sumMomentumXorYId (this->getStatistics().subscribeSum()),
      parameterPorosity(parameterPorosity_),
      flowDirection(flowDirection_)
{ }

template<typename T, template<typename U> class Descriptor>
void DarcyVelocityFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& latticeGas, ScalarField2D<T>& geometry )
{
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                Cell<T,Descriptor>& cell=latticeGas.get(iX,iY);

                // Compute density and momentum on fluid cells only.
                plint materialIndex = util::roundToInt(geometry.get(iX,iY));
                if (materialIndex > parameterPorosity) {  // Fluid cell (BGK)
                //if (materialIndex < parameterPorosity) {  // Fluid cell (BGK)  per babak geoemtry
                    Array<T,Descriptor<T>::d> gasMomentum;
                    momentTemplates<T,Descriptor>::get_j(cell,gasMomentum);
                    statistics.gatherSum(sumMomentumXorYId, gasMomentum[flowDirection]);
                    statistics.gatherSum(sumDensityId, cell.computeDensity());
                }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
DarcyVelocityFunctional2D<T,Descriptor>* DarcyVelocityFunctional2D<T,Descriptor>::clone() const
{
    return new DarcyVelocityFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
T DarcyVelocityFunctional2D<T,Descriptor>::getSumDensity() const {
    return this->getStatistics().getSum(sumDensityId);
}

template<typename T, template<typename U> class Descriptor>
T DarcyVelocityFunctional2D<T,Descriptor>::getSumMomentumXorY() const {
    return this->getStatistics().getSum(sumMomentumXorYId);
}

template<typename T, template<typename U> class Descriptor>
void DarcyVelocityFunctional2D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;
}



template<typename T, template<typename U> class Descriptor>
T computeDarcyVelocity(MultiBlockLattice2D<T,Descriptor>& latticeGas, MultiScalarField2D<T>& geometry,
                       Box2D domain , plint parameterPorosity, plint flowDirection )
{
    DarcyVelocityFunctional2D<T,Descriptor> functional(parameterPorosity,flowDirection);
    applyProcessingFunctional(functional, domain, latticeGas, geometry);
    T sumDensity = functional.getSumDensity();
    T sumMomentumY = functional.getSumMomentumXorY();
    return sumMomentumY / sumDensity;
}


///

//////////////////////  Drag term 2D

template< typename T,template<typename U1> class Descriptor>
DragTerm2D<T,Descriptor>::
        DragTerm2D(T h_, T latticeNu_)
    :  h(h_), latticeNu(latticeNu_)
{ }



template< typename T,template<typename U1> class Descriptor>
void DragTerm2D<T,Descriptor>::
                    process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice)

{
    
    enum {
        forceOffset = Descriptor<T>::ExternalField::forceBeginsAt,
    };


    for (plint iX=domain.x0; iX<=domain.x1; ++iX) 
    {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) 
        {


            Array<T,Descriptor<T>::d> vel(0.,0.);
            Array<T,Descriptor<T>::d> Drag_component(0.,0.);
            lattice.get(iX,iY).computeVelocity(vel);                
            Drag_component[0] = -(T(12)*latticeNu/(h*h))*vel[0];
            Drag_component[1] = -(T(12)*latticeNu/(h*h))*vel[1];

            T *DragTerm = lattice.get(iX,iY).getExternal(forceOffset);
            Drag_component.to_cArray(DragTerm);
        }
            
    }
}
template< typename T,template<typename U1> class Descriptor>
DragTerm2D<T,Descriptor>*
    DragTerm2D<T,Descriptor>::clone() const
{
    return new DragTerm2D<T,Descriptor>(*this);
}

template< typename T, template<typename U> class Descriptor>
void DragTerm2D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    
    modified[0] = modif::staticVariables; // fluid
}  











/*
template<typename T, template<typename U> class Descriptor>
AverageMomentumFunctional2D<T,Descriptor>::AverageMomentumFunctional2D(plint parameterPorosity_)
    : averageMomentumId (this->getStatistics().subscribeAverage()),
      parameterPorosity(parameterPorosity_)
{ }

template<typename T, template<typename U> class Descriptor>
void AverageMomentumFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& geometry )
{
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                Cell<T,Descriptor>& cell=lattice.get(iX,iY);

                // Compute momentum on fluid cells only.
                plint materialIndex = util::roundToInt(geometry.get(iX,iY));
                if (materialIndex>parameterPorosity) {  // Fluid cell (BGK)
                    Array<T,Descriptor<T>::d> momentum;
                    momentTemplates<T,Descriptor>::get_j(cell,momentum);
                    statistics.gatherAverage(averageMomentumId, momentum[1]);
                    statistics.incrementStats();
                }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
AverageMomentumFunctional2D<T,Descriptor>* AverageMomentumFunctional2D<T,Descriptor>::clone() const
{
    return new AverageMomentumFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
T AverageMomentumFunctional2D<T,Descriptor>::getAverageMomentum() const {
    return this->getStatistics().getAverage(averageMomentumId);
}

template<typename T, template<typename U> class Descriptor>
T computeAverageMomentum(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& geometry,
                         plint yPos, plint parameterPorosity )
{
    AverageMomentumFunctional2D<T,Descriptor> functional(parameterPorosity);
    Box2D domain(0, lattice.getNx()-1, yPos, yPos);
    applyProcessingFunctional(functional, domain, lattice, geometry);
    T averageMomentum = functional.getAverageMomentum();
    return averageMomentum;
}


template<typename T, template<typename U> class Descriptor>
AverageVelocityFunctional2D<T,Descriptor>::AverageVelocityFunctional2D(plint parameterPorosity_)
    : averageVelocityId (this->getStatistics().subscribeAverage()),
      parameterPorosity(parameterPorosity_)
{ }

template<typename T, template<typename U> class Descriptor>
void AverageVelocityFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& geometry )
{
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                Cell<T,Descriptor>& cell=lattice.get(iX,iY);

                // Compute momentum on fluid cells only.
                plint materialIndex = util::roundToInt(geometry.get(iX,iY));
                if (materialIndex>parameterPorosity) {  // Fluid cell (BGK)
                    Array<T,Descriptor<T>::d> velocity;
                    momentTemplates<T,Descriptor>::compute_uLb(cell,velocity);
                    statistics.gatherAverage(averageVelocityId, velocity[1]);
                    statistics.incrementStats();
                }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
AverageVelocityFunctional2D<T,Descriptor>* AverageVelocityFunctional2D<T,Descriptor>::clone() const
{
    return new AverageVelocityFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
T AverageVelocityFunctional2D<T,Descriptor>::getAverageVelocity() const {
    return this->getStatistics().getAverage(averageVelocityId);
}



template<typename T, template<typename U> class Descriptor>
T computeAverageVelocity(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& geometry,
                         plint yPos, plint parameterPorosity )
{
    AverageVelocityFunctional2D<T,Descriptor> functional(parameterPorosity);
    Box2D domain(0, lattice.getNx()-1, yPos, yPos);
    applyProcessingFunctional(functional, domain, lattice, geometry);
    T averageMomentum = functional.getAverageVelocity();
    return averageMomentum;
}*/

#endif
