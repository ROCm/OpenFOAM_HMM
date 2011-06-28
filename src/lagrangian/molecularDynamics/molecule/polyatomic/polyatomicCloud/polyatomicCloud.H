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
    Foam::polyatomicCloud

Description

SourceFiles
    polyatomicCloudI.H
    polyatomicCloud.C


\*---------------------------------------------------------------------------*/

#ifndef polyatomicCloud_H
#define polyatomicCloud_H

#include "Cloud.H"
#include "polyatomic.H"
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
                       Class polyatomicCloud Declaration
\*---------------------------------------------------------------------------*/

class polyatomicCloud
:
    public Cloud<polyatomic>
{

private:

    // Private data

        const polyMesh& mesh_;

        const potential& pot_;

        List<DynamicList<polyatomic*> > cellOccupancy_;

        InteractionLists<polyatomic> il_;

        List<polyatomic::constantProperties> constPropList_;

        Random rndGen_;


    // Private Member Functions

        void buildConstProps();

        void setSiteSizesAndPositions();

        //- Determine which polyatomics are in which cells
        void buildCellOccupancy();

        void calculatePairForce();

        inline void evaluatePair
        (
            polyatomic& molI,
            polyatomic& molJ
        );

        inline bool evaluatePotentialLimit
        (
            polyatomic& molI,
            polyatomic& molJ
        ) const;

        void calculateTetherForce();

        void calculateExternalForce();

        void removeHighEnergyOverlaps();

        void initialiseMolecules
        (
            const IOdictionary& mdInitialiseDict
        );

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

        label nSites() const;

        inline vector equipartitionLinearVelocity
        (
            scalar temperature,
            scalar mass
        );

        inline vector equipartitionAngularMomentum
        (
            scalar temperature,
            const polyatomic::constantProperties& cP
        );

        //- Disallow default bitwise copy construct
        polyatomicCloud(const polyatomicCloud&);

        //- Disallow default bitwise assignment
        void operator=(const polyatomicCloud&);


public:

    // Constructors

        //- Construct given mesh and potential references
        polyatomicCloud
        (
            const polyMesh& mesh,
            const potential& pot,
            bool readFields = true
        );

        //- Construct given mesh, potential and mdInitialiseDict
        polyatomicCloud
        (
            const polyMesh& mesh,
            const potential& pot,
            const IOdictionary& mdInitialiseDict,
            bool readFields = true
        );


    // Member Functions

        //- Evolve the polyatomics (move, calculate forces, control state etc)
        void evolve();

        void calculateForce();

        void applyConstraintsAndThermostats
        (
            const scalar targetTemperature,
            const scalar measuredTemperature
        );


        // Access

            inline const polyMesh& mesh() const;

            inline const potential& pot() const;

            inline const List<DynamicList<polyatomic*> >& cellOccupancy() const;

            inline const InteractionLists<polyatomic>& il() const;

            inline const List<polyatomic::constantProperties> constProps() const;

            inline const polyatomic::constantProperties&
                constProps(label id) const;

            inline Random& rndGen();


    // Member Operators

        //- Write polyatomic sites in XYZ format
        void writeXYZ(const fileName& fName) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyatomicCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
