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

#ifndef SINGLE_SPECIES_HETEROGENEOUS_CHEMICAL_REACTIONS_PROCESSING_FUNCTIONAL_2D_H
#define SINGLE_SPECIES_HETEROGENEOUS_CHEMICAL_REACTIONS_PROCESSING_FUNCTIONAL_2D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "atomicBlock/reductiveDataProcessingFunctional2D.h"
#include "atomicBlock/atomicContainerBlock2D.h"
#include "core/plbDebug.h"
#include "core/util.h"
#include "core/blockStatistics.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "finiteDifference/fdStencils1D.h"
#include "atomicBlock/atomicBlock2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"

//#include "offLattice/triangleBoundary2D.h"
// #include "algorithm/functions.h"
#include <Eigen3/Core>

namespace plb {

template<typename T>
class getPH2D : public BoxProcessingFunctional2D_SS<T,T>
{
public :
   
    virtual void process(Box2D domain ,ScalarField2D<T>& H,ScalarField2D<T>& PH);
    virtual getPH2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;

};

template<typename T, template<class U> class Descriptor1>
class CalculateSourceSinkAndUpdateSolidFractionProcessor2D : public BoxProcessingFunctional2D
{
public :
    CalculateSourceSinkAndUpdateSolidFractionProcessor2D(T Ke5_, T alpha_, T M_, T imposedConcentration_);
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);

    virtual CalculateSourceSinkAndUpdateSolidFractionProcessor2D<T,Descriptor1>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private :
    T Ke5, alpha, M, imposedConcentration ;
};
    
    


   
template< typename T, template<typename U> class AdvDiffDescriptor, template<typename U> class FluidDescriptor >
class UpdateFlagProcess2D : public BoxProcessingFunctional2D
{
public :
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);

    virtual UpdateFlagProcess2D<T,AdvDiffDescriptor,FluidDescriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;

};
    
    
    
    
    
    
template<typename T, template<class U> class Descriptor1>
class CalculateSourceSinkPerPrimariesSpeciesProcessor2D : public BoxProcessingFunctional2D
{
public :
    CalculateSourceSinkPerPrimariesSpeciesProcessor2D(T M_, T omegaAD_, T imposedConcentration1_, bool nonActiveNodesBounceBack_);
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);

    virtual CalculateSourceSinkPerPrimariesSpeciesProcessor2D<T,Descriptor1>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private :
    T M, omegaAD, imposedConcentration1;
    bool nonActiveNodesBounceBack;
};
    
    


template< typename T, template<typename U> class FluidDescriptor,template<typename U> class AdvDiffDescriptor>
class ChangeDynamicDissolutionPrecipitationProcessor2D : public BoxProcessingFunctional2D
{
public :
    ChangeDynamicDissolutionPrecipitationProcessor2D(T omega_, Array<T,2> vectorForce_, T iniRho_);
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);

    virtual ChangeDynamicDissolutionPrecipitationProcessor2D<T,FluidDescriptor,AdvDiffDescriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private :
    T omega ; 
    Array<T,2> vectorForce;
    T iniRho;
};
    

template< typename T,template<typename U> class Descriptor>
class UpdateSolidFractionAndVelocityInExternalFields2D : public BoxProcessingFunctional2D {
public:
//     UpdateSolidFractionAndVelocityInExternalFields2D(Array<T,3> velWall_, T rhoIni_, T omega_)
//         : velWall(velWall_), rhoIni(rhoIni_), omega(omega_)
//     { }
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual UpdateSolidFractionAndVelocityInExternalFields2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
// private:
//     Array<T,3> velWall;
//     T    rhoIni, omega;
};





template< typename T, template<typename U> class FluidDescriptor,template<typename U> class AdvDiffDescriptor>
class PassiveScalarProcessor2D : public BoxProcessingFunctional2D
{
public :

    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);

    virtual PassiveScalarProcessor2D<T,FluidDescriptor,AdvDiffDescriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;

};


template<typename T, template<class U> class AdDescriptor>
class ChrisDiricheletProcessor2D : public BoxProcessingFunctional2D
{
public :
    ChrisDiricheletProcessor2D(T imposedC_, T omega_);
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);

    virtual ChrisDiricheletProcessor2D<T,AdDescriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private :
    T imposedC;
    T omega ; 
};


template< typename T, template<typename U> class FluidDescriptor,template<typename U> class AdvDiffDescriptor>
class CorrectionTermAdvectionDiffusionBGKProcessor2D : public BoxProcessingFunctional2D
{
public :
    CorrectionTermAdvectionDiffusionBGKProcessor2D(T omega_);
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);

    virtual CorrectionTermAdvectionDiffusionBGKProcessor2D<T,FluidDescriptor,AdvDiffDescriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private :
    T omega ; 
};



}  // namespace plb

#endif  // CHEMICAL_REACTIONS_PROCESSING_FUNCTIONAL_2D_H