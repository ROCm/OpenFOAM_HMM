/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::MassTransferPhaseSystem

Description
    Class for mass transfer between phases

SourceFiles
    MassTransferPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef MassTransferPhaseSystem_H
#define MassTransferPhaseSystem_H

#include "phaseSystem.H"
#include "HashPtrTable.H"
#include "interfaceCompositionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class MassTransferPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class MassTransferPhaseSystem
:
    public BasePhaseSystem
{
public:

    // Public typedefs

     typedef
            HashTable
            <
                autoPtr<interfaceCompositionModel>,
                phasePairKey,
                phasePairKey::hash
            >
            massTransferModelTable;

protected:

    // Protected typedefs

       typedef
            HashPtrTable
            <
                volScalarField,
                phasePairKey,
                phasePairKey::hash
            >
            dmdtTable;


    // Protected data

        //- Overall inter-phase mass transfer rates [Kg/s]
        dmdtTable dmdt_;

        //- Mass transfer models
        massTransferModelTable massTransferModels_;


    // Protected memebers

        //- Calculate L between phases
        tmp<volScalarField> calculateL
        (
            const volScalarField& dmdtNetki,
            const phasePairKey& keyik,
            const phasePairKey& keyki,
            const volScalarField& T
        ) const;


public:

    // Constructors

        //- Construct from fvMesh
        MassTransferPhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~MassTransferPhaseSystem() = default;


    // Member Functions

        //- Return total interfacial mass flow rate
        tmp<volScalarField> dmdt(const phasePairKey& key) const;


    // Mass transfer functions

        //- Return the heat transfer matrix and fill dmdt for phases
        virtual  tmp<fvScalarMatrix> heatTransfer(const volScalarField& T);


        //- Calculate mass transfer for species
        virtual void massSpeciesTransfer
        (
            const phaseModel& phase,
            volScalarField::Internal& Su,
            volScalarField::Internal& Sp,
            const word speciesName
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MassTransferPhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
