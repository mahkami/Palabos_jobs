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

#ifndef POROUS_MEDIA_H
#define POROUS_MEDIA_H

#include "palabos2D.h"
using namespace plb;

template<typename T, template<typename U> class Descriptor>
void constructPorousMedia(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& texture,  plint parameterPorosity);

template<typename T, template<typename U> class Descriptor>
void initializeFlagStatus(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& flagStatus);

template<typename T, template<typename U> class Descriptor>
void cleanFlagStatus(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& flagStatus, plint iniRho, Array<T,2> u, T omega);

template<typename T, template<typename U> class Descriptor>
void initializeFlagStatusTest(MultiBlockLattice2D<T,Descriptor>& species1, MultiScalarField2D<T>& solidFraction, MultiScalarField2D<T>& flagStatus);

template<typename T, template<typename U> class Descriptor>
void initializeSolidFraction(MultiBlockLattice2D<T,Descriptor>& lattice, MultiScalarField2D<T>& solidFraction);

template<typename T, template<typename U> class Descriptor>
void perPorousMediaVisulatization(MultiScalarField2D<T>& texture,  MultiScalarField2D<T>& porousMedium,  plint parameterPorosity);

template<typename T, template<typename U> class Descriptor>
void initializeDynamicPrincipleSpecies(MultiBlockLattice2D<T,Descriptor>& species3, MultiBlockLattice2D<T,Descriptor>& species2,
                                       MultiBlockLattice2D<T,Descriptor>& species1, MultiScalarField2D<T>& porousMediumFlag, 
                                       bool nonActiveNodesBounceBack_, T iniConcentration1_, T iniConcentration2_, T iniConcentration3_  );

                                       
template<typename T, template<typename U> class Descriptor>
void initializeDynamicSinglepecies(MultiBlockLattice2D<T,Descriptor>& species1, MultiScalarField2D<T>& porousMediumFlag, 
                                   bool nonActiveNodesBounceBack, T imposedConcentration, T omegaAD);
                                       
template<typename T>
T computePorosity(MultiScalarField2D<T>& geometry, plint parameterPorosity, Box2D domain);
                       
template<typename T, template<typename U> class Descriptor>
T computePorosityFromFluidLattice(MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain);

template<typename T, template<typename U> class Descriptor>
T computeDarcyVelocity(MultiBlockLattice2D<T,Descriptor>& latticeGas, MultiScalarField2D<T>& geometry,
                       Box2D domain , plint parameterPorosity, plint flowDirection );

template<typename T, template<typename U> class Descriptor>
void porousMediaEvolution(MultiBlockLattice2D<T,Descriptor>& latticeGa, MultiScalarField2D<T>& porousMediumFlag);
// 
// template<typename T, template<typename U> class Descriptor>
// T computeAverageMomentum(MultiBlockLattice2D<T,Descriptor>& latticeGas, MultiScalarField2D<T>& geometry,
//                          plint yPos, plint parameterPorosity );





template<typename T, template<typename U> class Descriptor>
class ConstructPorousMediaFunctional2D : public BoxProcessingFunctional2D_LS<T,Descriptor,T>
{
public:
    ConstructPorousMediaFunctional2D(plint parameterPorosity_ );
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                            ScalarField2D<T>& geometry);
    virtual ConstructPorousMediaFunctional2D<T,Descriptor>* clone() const;
    // Envelope must be included, because a new dynamics is assigned to cells
    // (and must be assigned to envelope too).
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    plint parameterPorosity;
};



template<typename T, template<typename U> class Descriptor>
class CleaningFlagStatusFunctional2D : public BoxProcessingFunctional2D_LS<T,Descriptor,T>
{
public:
    CleaningFlagStatusFunctional2D(T iniRho_, Array<T,2> u_, T omega_ );
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,ScalarField2D<T>& flagStatus);
    virtual CleaningFlagStatusFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    T iniRho; 
    Array<T,2> u; 
    T omega;
};



template<typename T, template<typename U> class Descriptor>
class InitializeFlagStatusFunctional2D : public BoxProcessingFunctional2D_LS<T,Descriptor,T>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,ScalarField2D<T>& flagStatus);
    virtual InitializeFlagStatusFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
};



template<typename T, template<typename U> class Descriptor>
class InitializeFlagStatusFunctionalTest2D : public BoxProcessingFunctional2D
{
public:
    void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual InitializeFlagStatusFunctionalTest2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
};


template<typename T, template<typename U> class Descriptor>
class PorousMediaEvolution2D : public BoxProcessingFunctional2D_LS<T,Descriptor,T>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,ScalarField2D<T>& porousMediumFlag);
    virtual PorousMediaEvolution2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
};

template<typename T, template<typename U> class Descriptor>
class InitializeSolifFractionFunctional2D : public BoxProcessingFunctional2D_LS<T,Descriptor,T>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,ScalarField2D<T>& solidFraction);
    virtual InitializeSolifFractionFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
};

template<typename T>
class VisualizePorousMedium2D : public BoxProcessingFunctional2D_SS<T,T>
{
public:
    VisualizePorousMedium2D(plint parameterPorosity_ );
    virtual void process(Box2D domain, ScalarField2D<T>& texture,ScalarField2D<T>& porousMedium);
    virtual VisualizePorousMedium2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    plint parameterPorosity;
};


template<typename T, template<typename U> class Descriptor>
class InitializeDynamicsSpecies2D : public BoxProcessingFunctional2D
{
public:
    InitializeDynamicsSpecies2D(bool nonActiveNodesBounceBack_, T iniConcentration1_, T iniConcentration2_, T iniConcentration3_ );
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual InitializeDynamicsSpecies2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    bool nonActiveNodesBounceBack;
    T iniConcentration1, iniConcentration2, iniConcentration3;
};


template<typename T, template<typename U> class Descriptor>
class InitializeDynamicsSingleSpecies2D : public BoxProcessingFunctional2D
{
public:
    InitializeDynamicsSingleSpecies2D(bool nonActiveNodesBounceBack_, T imposedConcentration_, T omegaAD_ );
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> atomicBlocks);
    virtual InitializeDynamicsSingleSpecies2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    bool nonActiveNodesBounceBack;
    T imposedConcentration;
    T omegaAD;
};



template<typename T>
class ComputePorosityFunctional2D : public ReductiveBoxProcessingFunctional2D_S<T>
{
public:
    ComputePorosityFunctional2D(T parameterPorosity_);
    virtual void process(Box2D domain, ScalarField2D<T>& geometry);
    virtual ComputePorosityFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    
    plint getNumFluidCells() const;
private:
    plint numFluidCellId;
    T parameterPorosity;
};



template<typename T, template<typename U> class Descriptor>
class ComputePorosityFromFluidFunctional2D : public ReductiveBoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    ComputePorosityFromFluidFunctional2D(bool flag);
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual ComputePorosityFromFluidFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    
    plint getNumFluidCells() const;
private:
    plint numFluidCellId;
    bool flag;
};



template<typename T, template<typename U> class Descriptor>
class DarcyVelocityFunctional2D : public ReductiveBoxProcessingFunctional2D_LS<T,Descriptor,T>
{
public:
    DarcyVelocityFunctional2D(plint parameterPorosity_, plint flowDirection_);
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& latticeGas, ScalarField2D<T>& geometry);
    virtual DarcyVelocityFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    
    T getSumDensity() const;
    T getSumMomentumXorY() const;
private:
    plint sumDensityId;
    plint sumMomentumXorYId;
    plint parameterPorosity;
    plint flowDirection;
};










template<typename T>
class ComputeSigmaValueFunctional2D : public ReductiveBoxProcessingFunctional2D_S<T>
{
public:
    ComputeSigmaValueFunctional2D(plint flagValue_);
    virtual void process(Box2D domain, ScalarField2D<T>& flagStatus);
    virtual ComputeSigmaValueFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    
    plint getNumInterfaceCells() const;
private:
    plint numInterfaceCellId;
    plint flagValue;
};




// Calculate drag term to 2D approximation as proposed in literature
template< typename T, template<typename U> class Descriptor>
class DragTerm2D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public :
    DragTerm2D(T h_, T latticeNu_);
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual DragTerm2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private :
    T h ;
    T latticeNu; 
};





/*
template<typename T, template<typename U> class Descriptor>
class AverageMomentumFunctional2D : public ReductiveBoxProcessingFunctional2D_LS<T,Descriptor,T>
{
public:
    AverageMomentumFunctional2D(plint parameterPorosity_);
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& latticeGas, ScalarField2D<T>& geometry);
    virtual AverageMomentumFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    T getAverageMomentum() const;
private:
    plint averageMomentumId;
    plint parameterPorosity;
};

template<typename T, template<typename U> class Descriptor>
class AverageVelocityFunctional2D : public ReductiveBoxProcessingFunctional2D_LS<T,Descriptor,T>
{
public:
    AverageVelocityFunctional2D(plint parameterPorosity_);
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& latticeGas, ScalarField2D<T>& geometry);
    virtual AverageVelocityFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    T getAverageVelocity() const;
private:
    plint averageVelocityId;
    plint parameterPorosity;
};*/

#endif
