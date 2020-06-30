/* This file is part of the Palabos library.
 * Copyright (C) 2009, 2010 Jonas Latt
 * E-mail contact: jonas@lbmethod.org
 * The most recent release of Palabos can be downloaded at 
 * <http://www.lbmethod.org/palabos/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/** \file
 * Helper functions for domain initialization -- header file.
 */
#ifndef SINGLE_SPECIES_HETEROGENEOUS_CHEMICAL_REACTIONS_PROCESSING_FUNCTIONAL_2D_HH
#define SINGLE_SPECIES_HETEROGENEOUS_CHEMICAL_REACTIONS_PROCESSING_FUNCTIONAL_2D_HH

#include "singleSpeciesHeterogeneousChemicalReactionsProcessingFunctional2D.h"
// #include "chrisDiricheletBoundaryDynamics.h"
// #include "chrisDiricheletBoundaryDynamics.hh"
#include "core/plbDebug.h"
#include "core/util.h"
#include "core/blockStatistics.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "finiteDifference/fdStencils1D.h"
#include "atomicBlock/atomicBlock2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include <cmath>
#include <limits>
#include <Eigen3/Core>
#include <Eigen3/LU>

namespace plb {
    
///////////////  Erosion Process   ////////////////////////////////////




template<typename T>
void getPH2D<T>::process (
        Box2D domain, ScalarField2D<T>& H, ScalarField2D<T>& PH )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                PH.get(iX,iY) = -log10(H.get(iX,iY));
        }
    }
}

template<typename T>
getPH2D<T>* getPH2D<T>::clone() const {
    return new getPH2D<T>(*this);
}

template<typename T>
void getPH2D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::dataStructure;
}

///////////////////////////////////////////////////////////////////////////////////////


template<typename T, template<class U> class Descriptor1>
CalculateSourceSinkAndUpdateSolidFractionProcessor2D<T,Descriptor1>::
    CalculateSourceSinkAndUpdateSolidFractionProcessor2D(T Ke5_, T alpha_, T M_, T imposedConcentration_ ) :
        Ke5(Ke5_), alpha(alpha_), M(M_), imposedConcentration(imposedConcentration_)
{  }

template<typename T, template<class U> class Descriptor1>
void CalculateSourceSinkAndUpdateSolidFractionProcessor2D<T,Descriptor1>::
        processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks)
{
     
    BlockLattice2D<T,Descriptor1>& species1 =
            *dynamic_cast<BlockLattice2D<T,Descriptor1> *>(atomicBlocks[0]);
   
    ScalarField2D<T>& solidFraction =
            *dynamic_cast<ScalarField2D<T> *>(atomicBlocks[1]);
    
    ScalarField2D<T>& saturation_index =
            *dynamic_cast<ScalarField2D<T> *>(atomicBlocks[2]);
   
    ScalarField2D<T>& changeInSolidFraction =
            *dynamic_cast<ScalarField2D<T> *>(atomicBlocks[3]);
            
    ScalarField2D<T>& flagMedia =
            *dynamic_cast<ScalarField2D<T> *>(atomicBlocks[4]);
            
    Dot2D offset = computeRelativeDisplacement(species1, saturation_index);
    Dot2D offsetFM = computeRelativeDisplacement(species1, flagMedia);
    Dot2D offsetSF = computeRelativeDisplacement(species1, solidFraction);
    Dot2D offsetCSF = computeRelativeDisplacement(species1, changeInSolidFraction);
    
    
    enum {
        mySourceSink   = Descriptor1<T>::ExternalField::scalarBeginsAt,
    };

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint oX = iX + offset.x;
        plint fmX = iX + offsetFM.x;
        plint sfX = iX + offsetSF.x;
        plint csfX = iX + offsetCSF.x;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint oY = iY + offset.y;
            plint fmY = iY + offsetFM.y;
            plint sfY = iY + offsetSF.y;
            plint csfY = iY + offsetCSF.y;
            
            
             if (flagMedia.get(fmX,fmY) == (T)1)
             {
                 
                T localConcentration = species1.get(iX,iY).computeDensity();
//                 if (localConcentration > T(1))  // good working stuff ... but let's give it a try
//                      localConcentration = T(1);
//                 else if (localConcentration < T(0))
//                     localConcentration = T(0);
                
                
                 if (localConcentration < T(0))
                    localConcentration = T(0);
                
                saturation_index.get(oX,oY) = localConcentration/imposedConcentration; // define the Keq5
                T liquidFraction = (T)1-solidFraction.get(sfX,sfY);
                
                T newLiquidFraction = liquidFraction - alpha*M*(saturation_index.get(oX,oY) - (T)1);
                
                //T newLiquidFraction = liquidFraction - Ke5*M*imposedConcentration*(saturation_index.get(oX,oY) - (T)1);
                

                
               /* if( newLiquidFraction > 1. )
                {    
                    newLiquidFraction = 1.;
                }
                else */
                if ( newLiquidFraction < 0. )
                {
                    newLiquidFraction = 0.;
                }
                              
                T changeInLiquidFraction = (newLiquidFraction - liquidFraction);
                
                changeInSolidFraction.get(csfX,csfY) = -changeInLiquidFraction;
                
                
                if( newLiquidFraction > 1. )
                {    
                    newLiquidFraction = 1.;
                }
                
                solidFraction.get(sfX,sfY) = (T)1 - newLiquidFraction;
            }
            

        }
    }
}

template<typename T, template<class U> class Descriptor1>
CalculateSourceSinkAndUpdateSolidFractionProcessor2D<T,Descriptor1>*
    CalculateSourceSinkAndUpdateSolidFractionProcessor2D<T,Descriptor1>::clone() const
{
    return new CalculateSourceSinkAndUpdateSolidFractionProcessor2D<T,Descriptor1>(*this);
}

template<typename T, template<class U> class Descriptor1>
void CalculateSourceSinkAndUpdateSolidFractionProcessor2D<T,Descriptor1>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
     // IF Full dynamics object must be recreated, because this data processor
    //   re-attributes a new dynamics and acts on the bulk only. dynamicVariables
    
    modified[0] = modif::nothing; // species3
    modified[1] = modif::staticVariables; // solidFraction
    modified[2] = modif::staticVariables; // saturation_index
    modified[3] = modif::staticVariables; // change in solid Fraction
    modified[4] = modif::nothing; // flagMedia
}


// well, I guess that once I have the new solidFraction I can calculate the new flagStatus
// and then, as a function of the new Flag Status, play with the dynamics of the advectionDiffusion Schemes

// 1. calculate the sink-source 
// and 2. if needed re-update the the dynamic of the node


//=====================================================================================
//============== UpdateFlagProcess2D ===============
//=====================================================================================


// After each iteration for fluid flow, the flagStatus for the Reactions nodes is checked and reupdated
// flagStatus.get(iX,iY) = 0 means that the lattice node is liquid --> the reaction will be active then
// flagStatus.get(iX,iY) = 1 means that the lattice node is solid but its neighs are liquid --> the reaction will be active then
// flagStatus.get(iX,iY) = 2 means that the lattice node is solid and surrounded by solid nodes --> the reaction won't be active. 
//                         So either BounceBack or DiricheletPhaseFieldMethod is working on the node.
// Moreover, we keep track of the old status of the flagMedia. Fun for output purposes and helpful for reupdating the dynamic of the 
// the reaction nodes. flagMedia and flagMediaOld are used in  ChangeDynamicDissolutionPrecipitationProcessor2D dataProcessor.

template< typename T, template<typename U> class AdvDiffDescriptor, template<typename U> class FluidDescriptor >
void UpdateFlagProcess2D<T,AdvDiffDescriptor,FluidDescriptor>::
        processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks)
{
    
    
    BlockLattice2D<T,AdvDiffDescriptor>& species1 =
            *dynamic_cast<BlockLattice2D<T,AdvDiffDescriptor> *>(atomicBlocks[0]);
                           
    ScalarField2D<T>& solidFraction =
            *dynamic_cast<ScalarField2D<T> *>(atomicBlocks[1]);
           
    ScalarField2D<T>& changeInSolidFraction =
            *dynamic_cast<ScalarField2D<T> *>(atomicBlocks[2]);
            
    ScalarField2D<T>& flagMedia =
            *dynamic_cast<ScalarField2D<T> *>(atomicBlocks[3]);
            
    ScalarField2D<T>& flagMediaOld =
            *dynamic_cast<ScalarField2D<T> *>(atomicBlocks[4]);
    
    BlockLattice2D<T,FluidDescriptor>& lattice =
            *dynamic_cast<BlockLattice2D<T,FluidDescriptor> *>(atomicBlocks[5]);
            
    Dot2D relativeOffsetSF = computeRelativeDisplacement(species1,solidFraction);
    Dot2D relativeOffsetFM = computeRelativeDisplacement(species1,flagMedia);
    Dot2D relativeOffset = computeRelativeDisplacement(species1,lattice);
//     BounceBack<T,FluidDescriptor> BBdynamics;

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) 
    {
        //plint adX = relativeOffset.x + iX;
        
        
        plint sfX       = relativeOffsetSF.x    + iX;
        plint fmX       = relativeOffsetFM.x    + iX;
        
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) 
        {
            plint sfY   = relativeOffsetSF.y    + iY;
            plint fmY   = relativeOffsetFM.y    + iY;
            
            T liquidFraction = (T)1 - solidFraction.get(sfX,sfY);
            
            
            // from Chris code !!!!  we use the liquid fraction values in order to define 
            // the activated/disactivated state of the cellAD
            
            
            
            plint test = 0;  // normalAdvectionDiffusion BGK or AdvectionDiffusion with Source but zero sink/source
            
            if (liquidFraction <= 0)
                test = 2;
            else
            {
                plint tbb = 0; 
                for(plint iPop=1; iPop<AdvDiffDescriptor<T>::q; iPop++)
                {
                    plint iXneighSF = sfX+AdvDiffDescriptor<T>::c[iPop][0];
                    plint iYneighSF = sfY+AdvDiffDescriptor<T>::c[iPop][1];
                    
                    T neighLiquidFraction = (T)1-solidFraction.get(iXneighSF,iYneighSF);
                    
                    //if( liquidFraction < (T)1 && neighLiquidFraction == (T)1)
                    if( liquidFraction < (T)1 && liquidFraction > (T)0 && neighLiquidFraction == (T)1)
                        test=1; // Source Sink Active
                    if( liquidFraction > (T)0 && neighLiquidFraction == (T)0 )
                        test=1; // Source Sink active
    //                 if (liquidFraction < 0.5 && neighLiquidFraction < 0.5) //?? // which one better? I think the second one
    //                 if (liquidFraction < 0.5 && neighLiquidFraction < 1.) //?? THIS IS THE BEST ONE
    //                     tbb += 1;
                    
//                     if (liquidFraction <= 0. && neighLiquidFraction < 1.) //?? THIS IS THE BEST ONE
//                         tbb += 1;
                    
                }
            }
        
           
//             if(tbb == AdvDiffDescriptor<T>::q - 1 ) // mines rest velocity 
//                 test=2;     // this is for ChrisDiricheletBoundaryBGKdynamics
            
            
            // well. This is the end of Chris' test for the determination of the flagStatus
            // I update the status of the flagMedia field
            
            flagMediaOld.get(fmX,fmY) = flagMedia.get(fmX,fmY); // in order to keep track of the previous state of the lattice node;
            flagMedia.get(fmX,fmY) = (T)test; // 0 stands for liquid Node;
           
 
        
        }
    }
}

template< typename T, template<typename U> class AdvDiffDescriptor, template<typename U> class FluidDescriptor >
UpdateFlagProcess2D<T,AdvDiffDescriptor,FluidDescriptor>*
    UpdateFlagProcess2D<T,AdvDiffDescriptor,FluidDescriptor>::clone() const
{
    return new UpdateFlagProcess2D<T,AdvDiffDescriptor,FluidDescriptor>(*this);
}


template< typename T, template<typename U> class AdvDiffDescriptor, template<typename U> class FluidDescriptor >
void UpdateFlagProcess2D<T,AdvDiffDescriptor,FluidDescriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
     // IF Full dynamics object must be recreated, because this data processor
    //   re-attributes a new dynamics and acts on the bulk only. dynamicVariables
    
   
    modified[0] = modif::nothing; // species1
    modified[1] = modif::nothing; // solidFraction
    modified[2] = modif::nothing; // changeInSolidFraction
    modified[3] = modif::staticVariables; // flagMedia
    modified[4] = modif::staticVariables; // flagMediaOld
    modified[5] = modif::nothing; // flagMediaOld
}






template<typename T, template<class U> class Descriptor1>
CalculateSourceSinkPerPrimariesSpeciesProcessor2D<T,Descriptor1>::
    CalculateSourceSinkPerPrimariesSpeciesProcessor2D(T M_, T omegaAD_,T imposedConcentration1_, bool nonActiveNodesBounceBack_  ) :
        M(M_),omegaAD(omegaAD_), imposedConcentration1(imposedConcentration1_), nonActiveNodesBounceBack(nonActiveNodesBounceBack_)
{  }

template<typename T, template<class U> class Descriptor1>
void CalculateSourceSinkPerPrimariesSpeciesProcessor2D<T,Descriptor1>::
        processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks)
{
     
    BlockLattice2D<T,Descriptor1>& species1 =
            *dynamic_cast<BlockLattice2D<T,Descriptor1> *>(atomicBlocks[0]);
    
    ScalarField2D<T>& changeInSolidFraction =
            *dynamic_cast<ScalarField2D<T> *>(atomicBlocks[1]);
            
    ScalarField2D<T>& flagMedia =
            *dynamic_cast<ScalarField2D<T> *>(atomicBlocks[2]);
            
    ScalarField2D<T>& flagMediaOld =
            *dynamic_cast<ScalarField2D<T> *>(atomicBlocks[3]);
            
    Dot2D offset = computeRelativeDisplacement(species1, changeInSolidFraction);
    Dot2D offsetFM = computeRelativeDisplacement(species1, flagMedia);
    Dot2D offsetFMO = computeRelativeDisplacement(species1, flagMediaOld);
    
    enum {
        mySourceSink   = Descriptor1<T>::ExternalField::scalarBeginsAt,
        velOffset   = Descriptor1<T>::ExternalField::velocityBeginsAt,
    };
    
    Array<T,Descriptor1<T>::d> vel (0.,0.); // it must be a solid boundary. --> Is this true?
    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint oX = iX + offset.x;
        plint fmX = iX + offsetFM.x;
        plint fmoX = iX + offsetFMO.x;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint oY = iY + offset.y;
            plint fmY = iY + offsetFM.y;
            plint fmoY = iY + offsetFMO.y;
            //  ATTENZIONE !!!!!!!!!!!!!!!!!!!!!!!!!   OLD STUFFFFFFFFFFFFFFFFFFFF
            
//             
            ///////////////////////////  ATTENZIONE NEW STUFFFFFFFFFFFFFFFFFFF  ////////////
            if (flagMedia.get(fmX,fmY) == (T)1 ){
                *(species1.get(iX,iY).getExternal(mySourceSink)) = -changeInSolidFraction.get(oX,oY)/M;
            }
            else if(flagMedia.get(fmX,fmY) == (T)0 ){
                *(species1.get(iX,iY).getExternal(mySourceSink)) = (T)0;
            }
            else if(flagMedia.get(fmX,fmY) == (T)2 ){
                *(species1.get(iX,iY).getExternal(mySourceSink)) = (T)0;
            }
           
            
            if (flagMedia.get(fmX,fmY) != flagMediaOld.get(fmoX,fmoY))
            {    // we might have to change the Dynamic
                if (flagMedia.get(fmX,fmY) == (T)2 ) // either BounceBack or ChrisBoundaryDiricheletCondition test=2
                {   
                    if (nonActiveNodesBounceBack)
                    {
                        defineDynamics(species1, iX,iY, 
                                   new BounceBack<T,Descriptor1>(imposedConcentration1)); 
                    }
                    else
                    {   
                        defineDynamics(species1, iX,iY, 
                                    new NoDynamics<T,Descriptor1>(imposedConcentration1)); 
                                    
                    }
                }
                else if (flagMedia.get(fmX,fmY) == (T)1 ) // we become an active node on the boundary with sink/source active term
                {                                         // but please.. we are still solid 
                    if (flagMediaOld.get(fmoX,fmoY) == T(2)) // ANDREA NEW ... otherwise it is already AdvectionDiffusion (change in dynamics is only needed)
                                                             // when changing from bounceback to Advect diffusion --> dissolution.
                    {
                        defineDynamics(species1, iX,iY, new AdvectionDiffusionWithSourceBGKdynamics<T,Descriptor1>(omegaAD));  
                        iniCellAtEquilibrium(species1.get(iX,iY), imposedConcentration1, vel);
                    }
                    *(species1.get(iX,iY).getExternal(mySourceSink)) = -changeInSolidFraction.get(oX,oY)/M; // however this part is always correct
                    T *u1 = species1.get(iX,iY).getExternal(velOffset);
                    vel.to_cArray(u1);
                    
                    
                }    
                else if (flagMedia.get(fmX,fmY) == (T)0 ) // we are in the fluid
                {
                    // well in this case the SourceSinkTerm should be equal to zero ... So I keep the same 
                    // AdvectionDiffusionWithSourceBGKdynamics but I can check if the changeInSolidFraction is really = 0
                }
            }
            
            
              
        }
    }
}

template<typename T, template<class U> class Descriptor1>
CalculateSourceSinkPerPrimariesSpeciesProcessor2D<T,Descriptor1>*
    CalculateSourceSinkPerPrimariesSpeciesProcessor2D<T,Descriptor1>::clone() const
{
    return new CalculateSourceSinkPerPrimariesSpeciesProcessor2D<T,Descriptor1>(*this);
}

template<typename T, template<class U> class Descriptor1>
void CalculateSourceSinkPerPrimariesSpeciesProcessor2D<T,Descriptor1>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
     // IF Full dynamics object must be recreated, because this data processor
    //   re-attributes a new dynamics and acts on the bulk only. dynamicVariables
    
    modified[0] = modif::dataStructure; // species1
    modified[1] = modif::nothing; // change in solid fraction
    modified[2] = modif::nothing; // flagMedia
    modified[3] = modif::nothing; // flagMediaOld

}



///////////////////////////////////////////////////////////////////////


//=====================================================================================
//==============  ChangeDynamicDissolutionPrecipitationProcessor2D ===============
//=====================================================================================

template< typename T, template<typename U> class FluidDescriptor, template<typename U> class AdvDiffDescriptor>
ChangeDynamicDissolutionPrecipitationProcessor2D<T,FluidDescriptor,AdvDiffDescriptor>::
        ChangeDynamicDissolutionPrecipitationProcessor2D( T omega_, Array<T,2> vectorForce_, T iniRho_)
    :  omega(omega_), vectorForce(vectorForce_), iniRho(iniRho_)
{ }


template< typename T, template<typename U> class FluidDescriptor,template<typename U> class AdvDiffDescriptor>
void ChangeDynamicDissolutionPrecipitationProcessor2D<T,FluidDescriptor,AdvDiffDescriptor>::
        processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks)
{
 
    //std::cout << "Do I enter ChangeDynamicDissolutionPrecipitationProcessor2D " << std::endl;
    enum {
        appliedForce    = FluidDescriptor<T>::ExternalField::forceBeginsAt,
    };
    
    BlockLattice2D<T,FluidDescriptor>& fluidLattice =
            *dynamic_cast<BlockLattice2D<T,FluidDescriptor> *>(atomicBlocks[0]);
    
    BlockLattice2D<T,AdvDiffDescriptor>& species1 =
            *dynamic_cast<BlockLattice2D<T,AdvDiffDescriptor> *>(atomicBlocks[1]);
            
    ScalarField2D<T>& solidFraction =
            *dynamic_cast<ScalarField2D<T> *>(atomicBlocks[2]);
    
    ScalarField2D<T>& flagMedia =
            *dynamic_cast<ScalarField2D<T> *>(atomicBlocks[3]);
            
    ScalarField2D<T>& flagMediaOld =
            *dynamic_cast<ScalarField2D<T> *>(atomicBlocks[4]);
    
    Dot2D relativeOffset = computeRelativeDisplacement(fluidLattice,species1);
    Dot2D relativeOffsetPM = computeRelativeDisplacement(fluidLattice,solidFraction);
    Dot2D relativeOffsetFM = computeRelativeDisplacement(fluidLattice,flagMedia);
    
    BounceBack<T,FluidDescriptor> BBdynamics;
    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) 
    {
        plint adX = relativeOffset.x + iX;
        plint slX = relativeOffsetPM.x + iX;
        plint fmX = relativeOffsetFM.x + iX;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) 
        {
            plint adY = relativeOffset.y + iY;
            plint slY = relativeOffsetPM.y + iY;
            plint fmY = relativeOffsetFM.y + iY;
            
            //Cell<T,AdvDiffDescriptor>& cellAD = adLattice.get(iX,iY);
            Cell<T,FluidDescriptor>&   cellNS = fluidLattice.get(iX,iY);
            
            // you get the sinkSource Term --> the new Solid Fraction        
            
            // we could do that also using the flagMedia and the flagOldMedia info
            // probably it would be even better because it would be consistent with what I do for 
            // the chemical Species
            
            if( solidFraction.get(slX,slY) > 0.5 && cellNS.getDynamics().getId() != BBdynamics.getId()  )
            {   
                defineDynamics(fluidLattice, iX,iY,  new BounceBack<T,FluidDescriptor>());
            }
            else if( solidFraction.get(slX,slY) <= 0.5 &&  cellNS.getDynamics().getId() == BBdynamics.getId() )
            {                          
                defineDynamics(fluidLattice,iX,iY, new GuoExternalForceBGKdynamics<T,FluidDescriptor>(omega));  
                //fluidLattice.attributeDynamics(iX,iY,new GuoExternalForceBGKdynamics<T,FluidDescriptor>(omega));

                Array<T,FluidDescriptor<T>::d> u(0.,0.);
            
                iniCellAtEquilibrium(cellNS, iniRho, u);

                vectorForce.to_cArray( cellNS.getExternal(appliedForce) );   
                
            }
            
            // new way
            
//             if (flagMedia.get(fmX,fmY) != flagMediaOld.get(fmX,fmY))
//             {    // we might have to change the Dynamic
//                 if (flagMedia.get(fmX,fmY) == (T)1 && flagMediaOld.get(fmX,fmY) == (T)0 ) 
//                 {
//                     // we become a bounceBack node
//                      defineDynamics(fluidLattice, iX,iY,  new BounceBack<T,FluidDescriptor>());
//                 }
//                 else if (flagMedia.get(fmX,fmY) == (T)0 && flagMediaOld.get(fmX,fmY) == (T)1 ) // ChrisBoundaryDiricheletCondition test=2
//                 {
//                     // we become a liquid node
//                     defineDynamics(fluidLattice,iX,iY, new GuoExternalForceBGKdynamics<T,FluidDescriptor>(omega));  
//                 //fluidLattice.attributeDynamics(iX,iY,new GuoExternalForceBGKdynamics<T,FluidDescriptor>(omega));
// 
//                     Array<T,FluidDescriptor<T>::d> u(0.,0.);
//                 
//                     iniCellAtEquilibrium(cellNS, iniRho, u);
// 
//                     vectorForce.to_cArray( cellNS.getExternal(appliedForce) );   
//                 }
// //                 else if (flagMedia.get(fmX,fmY) == (T)2 && flagMediaOld.get(fmX,fmY) == (T)1 ) 
// //                 {
// //                     // or viceversa --> no change in dynamic for the fluid it should be still bounceBack
// //                 }
// //                 else if (flagMedia.get(fmX,fmY) == (T)2 && flagMediaOld.get(fmX,fmY) == (T)0 ) 
// //                 {
// //                     pcout << "I don't like that " << std::endl;
// //                 }
// //                 else if (flagMedia.get(fmX,fmY) == (T)0 && flagMediaOld.get(fmX,fmY) == (T)2 ) 
// //                 {
// //                     pcout << "I don't like that neither " << std::endl;
// //                 }
//                 
//             }
//             
//             // if there is no change in flagMedia there is no change in dynamic !!!
//             
//             // ciao ciao
         }
    }
}

template< typename T, template<typename U> class FluidDescriptor,template<typename U> class AdvDiffDescriptor>
ChangeDynamicDissolutionPrecipitationProcessor2D<T,FluidDescriptor,AdvDiffDescriptor>*
    ChangeDynamicDissolutionPrecipitationProcessor2D<T,FluidDescriptor,AdvDiffDescriptor>::clone() const
{
    return new ChangeDynamicDissolutionPrecipitationProcessor2D<T,FluidDescriptor,AdvDiffDescriptor>(*this);
}


template< typename T, template<typename U> class FluidDescriptor,template<typename U> class AdvDiffDescriptor>
void ChangeDynamicDissolutionPrecipitationProcessor2D<T,FluidDescriptor,AdvDiffDescriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
     // IF Full dynamics object must be recreated, because this data processor
    //   re-attributes a new dynamics and acts on the bulk only. dynamicVariables
    
    modified[0] = modif::dataStructure; // fluid  // not 100 % sure Though .. maybe it should work on bulkAndEnvelope?
    modified[1] = modif::nothing; // species1
    modified[2] = modif::nothing; // solidFraction
    modified[3] = modif::nothing; // species1
    modified[4] = modif::nothing; // solidFraction
}



//////////////////


template<typename T, template<typename U> class Descriptor>
void UpdateSolidFractionAndVelocityInExternalFields2D<T,Descriptor>::processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks)

{
    
    PLB_ASSERT(atomicBlocks.size() == 2);

    
    BlockLattice2D<T,Descriptor>* lattice = dynamic_cast<BlockLattice2D<T,Descriptor> *>(atomicBlocks[0]);
    PLB_ASSERT(lattice);

    ScalarField2D<T>* solidFraction = dynamic_cast<ScalarField2D<T>*>(atomicBlocks[1]);
    PLB_ASSERT(solidFraction);
        
//     TensorField2D<T,2>* wallVelocity = dynamic_cast<TensorField2D<T,3>*>(atomicBlocks[2]); 
//     PLB_ASSERT(wallVelocity);
//     
    
    Dot2D of_SF = computeRelativeDisplacement(*lattice, *solidFraction);
//     Dot2D of_WV = computeRelativeDisplacement(*lattice, *wallVelocity);    
    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint iXsf      = iX + of_SF.x;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint iYsf      = iY + of_SF.y;

  
                Cell<T,Descriptor>&   cell = lattice->get(iX,iY);
                *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt) = solidFraction->get(iXsf,iYsf); // solid fraction and velocity wall
                
                T *wallVelocityFromExternal = cell.getExternal(Descriptor<T>::ExternalField::jBeginsAt);
                Array<T,Descriptor<T>::d> velWall;

                velWall[0] = T();
                velWall[1] = T();

                velWall.to_cArray(wallVelocityFromExternal);
                

        }
    }
}

template< typename T,template<typename U> class Descriptor>
UpdateSolidFractionAndVelocityInExternalFields2D<T,Descriptor>* UpdateSolidFractionAndVelocityInExternalFields2D<T,Descriptor>::clone() const {
    return new UpdateSolidFractionAndVelocityInExternalFields2D<T,Descriptor>(*this);
}


template<typename T, template<typename U> class Descriptor>
void UpdateSolidFractionAndVelocityInExternalFields2D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{

    modified[0] = modif::dataStructure;     // lattice
    modified[1] = modif::nothing;           // solidFraction
//     modified[2] = modif::nothing;           // velocityWall
 
}







//// coupling with Fluid


template< typename T,template<typename U1> class FluidDescriptor, template<typename U2> class AdDescriptor>
void PassiveScalarProcessor2D<T,FluidDescriptor,AdDescriptor>::
                    processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks)
{

    enum {
        velOffset   = AdDescriptor<T>::ExternalField::velocityBeginsAt,
//         forceOffset = FluidDescriptor<T>::ExternalField::forceBeginsAt,
//          mySourceSink   = AdDescriptor<T>::ExternalField::scalarBeginsAt,
    };

    BlockLattice2D<T,FluidDescriptor>& fluid =
            *dynamic_cast<BlockLattice2D<T,FluidDescriptor> *>(atomicBlocks[0]);
            
    BlockLattice2D<T,AdDescriptor>& species1 =
            *dynamic_cast<BlockLattice2D<T,AdDescriptor> *>(atomicBlocks[1]);
   

    Dot2D relativeOffset = computeRelativeDisplacement(fluid,species1);
    
    BounceBack<T,FluidDescriptor> BBdynamics;
    
    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX)
    {
        plint adX = relativeOffset.x + iX;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY)
        {
            plint adY = relativeOffset.y + iY;
             
            Cell<T,FluidDescriptor>&   cellNS = fluid.get(iX,iY);
        
            // Velocity coupling
            if(cellNS.getDynamics().getId() != BBdynamics.getId() )
            {
                T *u1 = species1.get(adX,adY).getExternal(velOffset);

            
                Array<T,AdDescriptor<T>::d> vel;
                fluid.get(iX,iY).computeVelocity(vel);
                vel.to_cArray(u1);
            }
            
            
            

        }
    }
    
    
}

template< typename T,template<typename U1> class FluidDescriptor, template<typename U2> class AdDescriptor>
PassiveScalarProcessor2D<T,FluidDescriptor,AdDescriptor>*
    PassiveScalarProcessor2D<T,FluidDescriptor,AdDescriptor>::clone() const
{
    return new PassiveScalarProcessor2D<T,FluidDescriptor,AdDescriptor>(*this);
}

template< typename T, template<typename U> class FluidDescriptor,template<typename U> class AdvDiffDescriptor>
void PassiveScalarProcessor2D<T,FluidDescriptor,AdvDiffDescriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    
    modified[0] = modif::nothing; // fluid
    modified[1] = modif::dataStructure; // species1

}  



////////////////  CHRIS DIRICHELET PROCESSOR

template<typename T, template<class U> class AdDescriptor>
ChrisDiricheletProcessor2D<T,AdDescriptor>::
        ChrisDiricheletProcessor2D( T imposedC_, T  omega_)
    :  imposedC(imposedC_), omega(omega_)
{ }

template<typename T, template<class U> class AdDescriptor>
void ChrisDiricheletProcessor2D<T,AdDescriptor>::
                    processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks)
{
 
    BlockLattice2D<T,AdDescriptor>& species1 =
            *dynamic_cast<BlockLattice2D<T,AdDescriptor> *>(atomicBlocks[0]);
   
    ScalarField2D<T>& flag =
            *dynamic_cast<ScalarField2D<T> *>(atomicBlocks[1]);

    Dot2D relativeOffset = computeRelativeDisplacement(species1,flag);
    
//    BounceBack<T,FluidDescriptor> BBdynamics;
    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX)
    {
        plint adX = relativeOffset.x + iX;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY)
        {
            plint adY = relativeOffset.y + iY;
            Cell<T,AdDescriptor>&   cell = species1.get(iX,iY); 
            
            if (flag.get(adX,adY) == 2 ) // you could also make sure that the Dynamic is NoDynamic
            {
                
                T rhoBar = T();
                for(plint iPop = 0; iPop<AdDescriptor<T>::q; iPop++)
                {
                    rhoBar += cell[iPop];
                }
                
                for (plint iPop=0; iPop < AdDescriptor<T>::q; ++iPop) {
                    cell[iPop] *= ((T)1-omega); 
                    cell[iPop] -= ((T)1-omega)*AdDescriptor<T>::t[iPop]*rhoBar; // feq*(1-omega)
                    cell[iPop] += AdDescriptor<T>::t[iPop]* (imposedC - 1.);
                }
                
            }
            
        }
    }
}

template<typename T, template<class U> class AdDescriptor>
ChrisDiricheletProcessor2D<T,AdDescriptor>*
    ChrisDiricheletProcessor2D<T,AdDescriptor>::clone() const
{
    return new ChrisDiricheletProcessor2D<T,AdDescriptor>(*this);
}

template<typename T, template<class U> class AdDescriptor>
void ChrisDiricheletProcessor2D<T,AdDescriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    
    modified[0] = modif::dynamicVariables; // fluid
    modified[1] = modif::nothing; // species1

}  

//////////////////////  correction term AD

template< typename T,template<typename U1> class FluidDescriptor, template<typename U2> class AdDescriptor>
CorrectionTermAdvectionDiffusionBGKProcessor2D<T,FluidDescriptor,AdDescriptor>::
        CorrectionTermAdvectionDiffusionBGKProcessor2D(T  omega_)
    :  omega(omega_)
{ }


template< typename T,template<typename U1> class FluidDescriptor, template<typename U2> class AdDescriptor>
void CorrectionTermAdvectionDiffusionBGKProcessor2D<T,FluidDescriptor,AdDescriptor>::
                    processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks)
{
    enum {
        velOffset               = AdDescriptor<T>::ExternalField::velocityBeginsAt,
        correctionTermAD        = AdDescriptor<T>::ExternalField::correctionTermBeginsAt,
//         forceOffset             = FluidDescriptor<T>::ExternalField::forceBeginsAt,
    };

    BlockLattice2D<T,FluidDescriptor>& fluid =
            *dynamic_cast<BlockLattice2D<T,FluidDescriptor> *>(atomicBlocks[0]);
            
    BlockLattice2D<T,AdDescriptor>& species1 =
            *dynamic_cast<BlockLattice2D<T,AdDescriptor> *>(atomicBlocks[1]);
            
         
    BounceBack<T,FluidDescriptor> BBdynamics;
    
    Dot2D relativeOffset = computeRelativeDisplacement(fluid,species1);

    
    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX)
    {
        plint adX = relativeOffset.x + iX;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY)
        {
            plint adY = relativeOffset.y + iY;

            // Let's take tje velocity at time t --> u(t) (same for the three Species)
            T *u1 = species1.get(adX,adY).getExternal(velOffset);
            
            // Let's take the velocity at time u(t+1) (same for the three species)
            Array<T,FluidDescriptor<T>::d> vel(0.,0.);
            if (fluid.get(iX,iY).getDynamics().getId() != BBdynamics.getId())
                fluid.get(iX,iY).computeVelocity(vel);
            
            
            ////////////////////  Let's calculate the correction term
            
            
            // geq0 change for different species.
            
            T fullRho_1species = species1.get(adX,adY).computeDensity();
            T rhoBar_1species  = fullRho_1species - 1.;
            T gADeq0_1species  =  AdDescriptor<T>::t[0]*rhoBar_1species;
            
            // non eqPart rest Particles. It changes for different species
            T nonEqRestParticleTerm_1stSpecies = species1.get(adX,adY)[0] - gADeq0_1species;
            nonEqRestParticleTerm_1stSpecies  *= omega;
            nonEqRestParticleTerm_1stSpecies  /= AdDescriptor<T>::t[0];
            
            Array<T,FluidDescriptor<T>::d> correctionTermAD_1stSpecies(0.,0.);
            for(plint i = 0; i < FluidDescriptor<T>::d; i++)
            {
                correctionTermAD_1stSpecies[i] = fullRho_1species*( (vel[i]-u1[i]) - nonEqRestParticleTerm_1stSpecies*vel[i] );
                correctionTermAD_1stSpecies[i] *= AdDescriptor<T>::invCs2*(T(1) - 0.5*omega); 
            }
           
            ////////////////////  End  calculation the correction term
            
            //////////////////  Store correction AD in externalField AD
            
            T *correctionTerm1 = species1.get(adX,adY).getExternal(correctionTermAD);
            correctionTermAD_1stSpecies.to_cArray(correctionTerm1);
            
         

        }
    }
}

template< typename T,template<typename U1> class FluidDescriptor, template<typename U2> class AdDescriptor>
CorrectionTermAdvectionDiffusionBGKProcessor2D<T,FluidDescriptor,AdDescriptor>*
    CorrectionTermAdvectionDiffusionBGKProcessor2D<T,FluidDescriptor,AdDescriptor>::clone() const
{
    return new CorrectionTermAdvectionDiffusionBGKProcessor2D<T,FluidDescriptor,AdDescriptor>(*this);
}

template< typename T, template<typename U> class FluidDescriptor,template<typename U> class AdvDiffDescriptor>
void CorrectionTermAdvectionDiffusionBGKProcessor2D<T,FluidDescriptor,AdvDiffDescriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    
    modified[0] = modif::nothing; // fluid
    modified[1] = modif::dataStructure; // species1
}  




}
#endif  //CHEMICAL_REACTIONS_PROCESSING_FUNCTIONAL_2D_HH