/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::initialPointsMethod

Description
    Abstract base class for generating initial points for a conformalVoronoiMesh

SourceFiles
    initialPointsMethod.C

\*---------------------------------------------------------------------------*/

#ifndef initialPointsMethod_H
#define initialPointsMethod_H

#include "point.H"
#include "conformalVoronoiMesh.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class initialPointsMethod Declaration
\*---------------------------------------------------------------------------*/

class initialPointsMethod
:
    public dictionary
{

protected:

    // Protected data

        //- Reference to the conformalVoronoiMesh holding this cvControls object
        const conformalVoronoiMesh& cvMesh_;

        //- Method details dictionary
        dictionary detailsDict_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        initialPointsMethod(const initialPointsMethod&);

        //- Disallow default bitwise assignment
        void operator=(const initialPointsMethod&);


public:

    //- Runtime type information
    TypeName("initialPointsMethod");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            initialPointsMethod,
            dictionary,
            (
                const dictionary& initialPointsDict,
                const conformalVoronoiMesh& cvMesh
            ),
            (initialPointsDict, cvMesh)
        );


    // Constructors

        //- Construct from components
        initialPointsMethod
        (
            const word& type,
            const dictionary& initialPointsDict,
            const conformalVoronoiMesh& cvMesh
        );


    // Selectors

        //- Return a reference to the selected initialPointsMethod
        static autoPtr<initialPointsMethod> New
        (
            const dictionary& initialPointsDict,
            const conformalVoronoiMesh& cvMesh
        );


    //- Destructor
    virtual ~initialPointsMethod();


    // Member Functions

        //- Const access to the details dictionary
        const dictionary& detailsDict() const
        {
            return detailsDict_;
        }

        //- Return the initial points for the conformalVoronoiMesh
        virtual std::vector<Vb::Point> initialPoints() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
