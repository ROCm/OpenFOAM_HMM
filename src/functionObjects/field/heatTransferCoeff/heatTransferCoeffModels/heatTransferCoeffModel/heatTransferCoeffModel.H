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

Namespace
    Foam::heatTransferCoeffModels

Description
    A namespace for various heat transfer coefficient model implementations.

Class
    Foam::heatTransferCoeffModel

Description
    An abstract base class for heat transfer coeffcient models.

SourceFiles
    heatTransferCoeffModel.C
    heatTransferCoeffModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransferCoeffModel_H
#define heatTransferCoeffModel_H

#include "dictionary.H"
#include "HashSet.H"
#include "volFields.H"

#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                   Class heatTransferCoeffModel Declaration
\*---------------------------------------------------------------------------*/

class heatTransferCoeffModel
{

    // Private Member Functions

        //- Disallow copy construct
        heatTransferCoeffModel(const heatTransferCoeffModel&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const heatTransferCoeffModel&) = delete;


protected:

    // Protected data

        const fvMesh& mesh_;

        const word TName_;

        labelHashSet patchSet_;

        word qrName_;

        tmp<FieldField<Field, scalar>> q() const;

        //- Set the heat transfer coefficient
        virtual void htc(volScalarField& htc) = 0;


public:

    //- Runtime type information
    TypeName("heatTransferCoeffModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            heatTransferCoeffModel,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh,
                const word& TName
            ),
            (dict, mesh, TName)
        );


    // Selectors

        //- Return a reference to the selected heat transfer coefficicent model
        static autoPtr<heatTransferCoeffModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& TName
        );


    // Constructors

        //- Construct from components
        heatTransferCoeffModel
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& TName
        );


    //- Destructor
    virtual ~heatTransferCoeffModel()
    {}


    // Member Functions

        //- Read from dictionary
        virtual bool read(const dictionary& dict);

        virtual bool calc(volScalarField& result);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
