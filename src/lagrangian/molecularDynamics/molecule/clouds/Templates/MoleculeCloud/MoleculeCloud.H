/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 OpenCFD Ltd.
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
    Foam::MoleculeCloud

Description

SourceFiles
    MoleculeCloudI.H
    MoleculeCloud.C


\*---------------------------------------------------------------------------*/

#ifndef MoleculeCloud_H
#define MoleculeCloud_H

#include "Cloud.H"
#include "moleculeCloud.H"
#include "IOdictionary.H"
#include "potential.H"
#include "InteractionLists.H"
#include "labelVector.H"
#include "Random.H"
#include "fileName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class MoleculeCloud Declaration
\*---------------------------------------------------------------------------*/

template<class MoleculeType>
class MoleculeCloud
:
    public Cloud<MoleculeType>,
    public moleculeCloud
{

private:

    // Private data

        //-
        const polyMesh& mesh_;

        //-
        const potential& pot_;

        //-
        List<DynamicList<MoleculeType*> > cellOccupancy_;

        //-
        InteractionLists<MoleculeType> il_;

        //-
        List<typename MoleculeType::constantProperties> constPropList_;

        //-
        Random rndGen_;


    // Private Member Functions

        //-
        void buildConstProps();

        //-
        void setSiteSizesAndPositions();

        //- Determine which molecules are in which cells
        void buildCellOccupancy();

        //-
        void calculatePairForce();

        //-
        inline void evaluatePair
        (
            MoleculeType& molI,
            MoleculeType& molJ
        );

        //-
        inline bool evaluatePotentialLimit
        (
            MoleculeType& molI,
            MoleculeType& molJ
        ) const;

        //-
        void calculateTetherForce();

        //-
        void calculateExternalForce();

        //-
        void removeHighEnergyOverlaps();

        //-
        void initialiseMolecules(const dictionary& mdInitialiseDict);

        //-
        void createMolecule
        (
            const point& position,
            label cell,
            label tetFace,
            label tetPt,
            label id,
            bool tethered,
            scalar temperature,
            const vector& bulkVelocity
        );

        //-
        label nSites() const;

        //- Disallow default bitwise copy construct
        MoleculeCloud(const MoleculeCloud&);

        //- Disallow default bitwise assignment
        void operator=(const MoleculeCloud&);


public:

    // Constructors

        //- Construct given mesh and potential references
        MoleculeCloud
        (
            const word& cloudName,
            const polyMesh& mesh,
            const potential& pot,
            bool readFields = true
        );

        //- Construct given mesh, potential and mdInitialiseDict
        MoleculeCloud
        (
            const word& cloudName,
            const polyMesh& mesh,
            const potential& pot,
            const dictionary& mdInitialiseDict,
            bool readFields = true
        );


    // Member Functions

        //- Evolve the molecules (move, calculate forces, control state etc)
        void evolve();

        //-
        void calculateForce();

        //- Print cloud information
        void info();


        // Access

            //-
            inline const polyMesh& mesh() const;

            //-
            inline const potential& pot() const;

            //-
            inline const List<DynamicList<MoleculeType*> >&
            cellOccupancy() const;

            //-
            inline const InteractionLists<MoleculeType>& il() const;

            //-
            inline const List<typename MoleculeType::constantProperties>
            constProps() const;

            //-
            inline const typename MoleculeType::constantProperties&
            constProps(label id) const;

            //-
            inline Random& rndGen();

            //-
            inline vector equipartitionLinearVelocity
            (
                scalar temperature,
                scalar mass
            );

            //-
            inline vector equipartitionAngularMomentum
            (
                scalar temperature,
                const typename MoleculeType::constantProperties& cP
            );


    // Member Operators

        //- Write molecule sites in XYZ format
        void writeXYZ(const fileName& fName) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MoleculeCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MoleculeCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
