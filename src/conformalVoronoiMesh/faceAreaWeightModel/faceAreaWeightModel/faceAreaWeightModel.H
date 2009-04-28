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
    Foam::faceAreaWeightModel

Description
    Abstract base class for providing faceAreaWeight values to the cell motion
    controller based on an argument faceAreaFraction value, typically the ratio
    of the area of the face in question to the area of a square face with a side
    length of the local target cell size.

SourceFiles
    faceAreaWeightModel.C

\*---------------------------------------------------------------------------*/

#ifndef faceAreaWeightModel_H
#define faceAreaWeightModel_H

#include "point.H"
#include "conformalVoronoiMesh.H"
#include "dictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class faceAreaWeightModel Declaration
\*---------------------------------------------------------------------------*/

class faceAreaWeightModel
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
        faceAreaWeightModel(const faceAreaWeightModel&);

        //- Disallow default bitwise assignment
        void operator=(const faceAreaWeightModel&);


public:

    //- Runtime type information
    TypeName("faceAreaWeightModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            faceAreaWeightModel,
            dictionary,
            (
                const dictionary& faceAreaWeightDict,
                const conformalVoronoiMesh& cvMesh
            ),
            (faceAreaWeightDict, cvMesh)
        );


    // Constructors

        //- Construct from components
        faceAreaWeightModel
        (
            const word& type,
            const dictionary& faceAreaWeightDict,
            const conformalVoronoiMesh& cvMesh
        );


    // Selectors

        //- Return a reference to the selected faceAreaWeightModel
        static autoPtr<faceAreaWeightModel> New
        (
            const dictionary& faceAreaWeightDict,
            const conformalVoronoiMesh& cvMesh
        );


    //- Destructor
    virtual ~faceAreaWeightModel();


    // Member Functions

        //- Const access to the coeffs dictionary
        const dictionary& coeffDict() const
        {
            return coeffDict_;
        }

        //- Return the current faceAreaWeight coefficient
        virtual scalar faceAreaWeight(scalar faceAreaFraction) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
