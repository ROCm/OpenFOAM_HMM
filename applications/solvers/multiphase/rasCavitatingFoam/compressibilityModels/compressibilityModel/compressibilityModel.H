/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    Foam::compressibilityModel

Description

SourceFiles
    compressibilityModel.C
    newCompressibilityModel.C

\*---------------------------------------------------------------------------*/

#ifndef compressibilityModel_H
#define compressibilityModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class compressibilityModel Declaration
\*---------------------------------------------------------------------------*/

class compressibilityModel
{

protected:

    // Protected data

        dictionary compressibilityProperties_;

        volScalarField psi_;
        const volScalarField& gamma_;


    // Private Member Functions

        //- Disallow copy construct
        compressibilityModel(const compressibilityModel&);

        //- Disallow default bitwise assignment
        void operator=(const compressibilityModel&);


public:

    //- Runtime type information
    TypeName("compressibilityModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            compressibilityModel,
            dictionary,
            (
                const dictionary& compressibilityProperties,
                const volScalarField& gamma
            ),
            (compressibilityProperties, gamma)
        );


    // Selectors

        //- Return a reference to the selected compressibility model
        static autoPtr<compressibilityModel> New
        (
            const dictionary& compressibilityProperties,
            const volScalarField& gamma
        );


    // Constructors

        //- Construct from components
        compressibilityModel
        (
            const dictionary& compressibilityProperties,
            const volScalarField& gamma
        );


    // Destructor

        virtual ~compressibilityModel()
        {}


    // Member Functions

        //- Return the phase transport properties dictionary
        const dictionary& compressibilityProperties() const
        {
            return compressibilityProperties_;
        }

        //- Return the compressibility
        const volScalarField& psi() const
        {
            return psi_;
        }

        //- Correct the compressibility
        virtual void correct() = 0;

        //- Read compressibilityProperties dictionary
        virtual bool read(const dictionary& compressibilityProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
