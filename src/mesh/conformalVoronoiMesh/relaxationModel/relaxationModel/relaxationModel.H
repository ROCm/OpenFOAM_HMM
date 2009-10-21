/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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
    Foam::relaxationModel

Description
    Abstract base class for providing relaxation values to the cell motion
    controller

SourceFiles
    relaxationModel.C

\*---------------------------------------------------------------------------*/

#ifndef relaxationModel_H
#define relaxationModel_H

#include "point.H"
#include "conformalVoronoiMesh.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class relaxationModel Declaration
\*---------------------------------------------------------------------------*/

class relaxationModel
:
    public dictionary
{

protected:

    // Protected data

        //- Reference to the conformalVoronoiMesh holding this cvControls object
        const conformalVoronoiMesh& cvMesh_;

        //- Method coeffs dictionary
        dictionary coeffDict_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        relaxationModel(const relaxationModel&);

        //- Disallow default bitwise assignment
        void operator=(const relaxationModel&);


public:

    //- Runtime type information
    TypeName("relaxationModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            relaxationModel,
            dictionary,
            (
                const dictionary& relaxationDict,
                const conformalVoronoiMesh& cvMesh
            ),
            (relaxationDict, cvMesh)
        );


    // Constructors

        //- Construct from components
        relaxationModel
        (
            const word& type,
            const dictionary& relaxationDict,
            const conformalVoronoiMesh& cvMesh
        );


    // Selectors

        //- Return a reference to the selected relaxationModel
        static autoPtr<relaxationModel> New
        (
            const dictionary& relaxationDict,
            const conformalVoronoiMesh& cvMesh
        );


    //- Destructor
    virtual ~relaxationModel();


    // Member Functions

        //- Const access to the coeffs dictionary
        const dictionary& coeffDict() const
        {
            return coeffDict_;
        }

        //- Return the current relaxation coefficient
        virtual scalar relaxation() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
