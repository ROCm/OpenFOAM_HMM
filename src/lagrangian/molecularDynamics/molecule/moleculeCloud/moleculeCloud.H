/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::moleculeCloud

Description

SourceFiles
    moleculeCloud.C


\*---------------------------------------------------------------------------*/

#ifndef moleculeCloud_H
#define moleculeCloud_H

#include "Cloud.H"
#include "molecule.H"
#include "IOdictionary.H"
#include "pairPotentialList.H"
#include "tetherPotentialList.H"
#include "interactionLists.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class moleculeCloud Declaration
\*---------------------------------------------------------------------------*/

class moleculeCloud
:
    public Cloud<molecule>
{

private:

    // Private data

        const polyMesh& mesh_;

            // MD solution parameters

            List<DynamicList<molecule*> > cellOccupancy_;

            interactionLists il_;

            List<molecule::constantProperties> constPropList_;;

    // Move to potential class and store reference to it ->
            scalar potentialEnergyLimit_;

            labelList removalOrder_;

            pairPotentialList pairPotentials_;

            tetherPotentialList tetherPotentials_;

            vector gravity_;
    // <- Move to potential class and store reference to it

        // Private Member Functions

            void applyConstraintsAndThermostats
            (
                const scalar targetTemperature,
                const scalar measuredTemperature
            );

            void calculateForce();

            void calculatePairForce();

            void calculateTetherForce();

            void calculateExternalForce();

            void removeHighEnergyOverlaps();

            //- Disallow default bitwise copy construct
            moleculeCloud(const moleculeCloud&);

            //- Disallow default bitwise assignment
            void operator=(const moleculeCloud&);


public:

    // Static data members

        static scalar kb;

        static scalar elementaryCharge;

        static scalar vacuumPermittivity;

    // Constructors

        //- Construct given mesh
        moleculeCloud(const polyMesh&);

        //- Construct given polyMesh and fields of components.  Intended for use
        //- by the molConfig utility
        // moleculeCloud
        // (
        //     const polyMesh& mesh,
        //     label nMol,
        //     const labelField& id,
        //     const scalarField& mass,
        //     const vectorField& positions,
        //     const labelField& cells,
        //     const vectorField& U,
        //     const vectorField& A,
        //     const labelField& tethered,
        //     const vectorField& tetherPositions
        // );


    // Member Functions

    //- Determine which molecules are in which cells
        void buildCellOccupancy();

        //- Evolve the molecules (move, calculate forces, control state etc)
        void evolve();

        //- Integrate the equations of motion using algorithm selected at
        //  runtime from a dictionary. This will also call the function
        //  to calculate the intermolecular forces (calculatePairForce()).



        // Access

            inline const polyMesh& mesh() const;

            inline scalar potentialEnergyLimit() const;

            inline const labelList& removalOrder() const;

            inline const pairPotentialList& pairPotentials() const;

            inline const tetherPotentialList& tetherPotentials() const;

            inline const vector& gravity() const;

            inline const List<DynamicList<molecule*> >& cellOccupancy() const;



    // Member Operators

        //- Write fields
        void writeFields() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "moleculeCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
