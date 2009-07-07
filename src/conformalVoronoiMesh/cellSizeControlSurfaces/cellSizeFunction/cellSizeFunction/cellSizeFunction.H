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
    Foam::cellSizeFunction

Description
    Abstract base class for specifying target cell sizes

SourceFiles
    cellSizeFunction.C

\*---------------------------------------------------------------------------*/

#ifndef cellSizeFunction_H
#define cellSizeFunction_H

#include "point.H"
#include "conformalVoronoiMesh.H"
#include "searchableSurface.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class cellSizeFunction Declaration
\*---------------------------------------------------------------------------*/

class cellSizeFunction
:
    public dictionary
{

public:

    //- Runtime type information
    TypeName("cellSizeFunction");

    //- Surface side mode
    enum sideMode
    {
        INSIDE,         // Control inside the surface
        OUTSIDE,        // Control outside the surface
        BOTHSIDES       // Control on both sides of a surface
    };


protected:

    // Protected data

        //- Reference to the conformalVoronoiMesh holding this cvs object
        const conformalVoronoiMesh& cvMesh_;

        //- Reference to the searchableSurface that cellSizeFunction
        //  relates to
        const searchableSurface& surface_;

        //- Method details dictionary
        dictionary coeffsDict_;

        //- Mode of size specification, i.e. inside, outside or bothSides
        sideMode sideMode_;

        //- Priority of this cellSizeFunction
        label priority_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        cellSizeFunction(const cellSizeFunction&);

        //- Disallow default bitwise assignment
        void operator=(const cellSizeFunction&);


public:

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            cellSizeFunction,
            dictionary,
            (
                const dictionary& cellSizeFunctionDict,
                const conformalVoronoiMesh& cvMesh,
                const searchableSurface& surface
            ),
            (cellSizeFunctionDict, cvMesh, surface)
        );


    // Constructors

        //- Construct from components
        cellSizeFunction
        (
            const word& type,
            const dictionary& cellSizeFunctionDict,
            const conformalVoronoiMesh& cvMesh,
            const searchableSurface& surface
        );


    // Selectors

        //- Return a reference to the selected cellSizeFunction
        static autoPtr<cellSizeFunction> New
        (
            const dictionary& cellSizeFunctionDict,
            const conformalVoronoiMesh& cvMesh,
            const searchableSurface& surface
        );


    //- Destructor
    virtual ~cellSizeFunction();


    // Member Functions

        //- Const access to the details dictionary
        inline const dictionary& coeffsDict() const
        {
            return coeffsDict_;
        }

        inline label priority() const
        {
            return priority_;
        }

        //- Modify scalar argument to the cell size specified by function.
        //  Return a boolean specifying if the function was used, i.e. false if
        //  the point was not in range of the surface for a spatially varying
        //  size.
        virtual bool cellSize(const point& pt, scalar& size) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
