/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::wallBoundedStreamLine

Group
    grpFieldFunctionObjects

Description
    This function object generates streamline data by sampling a set of
    user-specified fields along a particle track, transported by a
    user-specified velocity field, constrained to a patch.

    Example of function object specification:
    \verbatim
    wallBoundedStreamLine1
    {
        type            wallBoundedStreamLine;
        functionObjectLibs ("libfieldFunctionObjects.so");
        ...
        setFormat       vtk;
        UName           UNear;
        trackForward    yes;
        fields
        (
            UNear
            p
        );
        lifeTime        10000;
        trackLength     1e-3;
        bounds          (0.2 -10 -10)(0.22 10 10);
        cloudName       particleTracks;
        seedSampleSet   patchSeed;
        patchSeedCoeffs
        {
            type        patchSeed;
            patches     (wall);
            axis        x;
            maxPoints   20000;
        }
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: wallBoundedStreamLine| yes |
        setFormat    | output data type        | yes         |
        UName        | tracking velocity field name | yes    |
        fields       | fields to sample        | yes         |
        lifetime     | maximum number of particle tracking steps | yes |
        trackLength  | tracking segment length | no          |
        cloudName    | cloud name to use       | yes         |
        log          | Log to standard output  | no          | yes
        bounds       | Bounding box to trim tracks | no | greatBox
        seedSampleSet| seeding method (see below)| yes       |
    \endtable

    \linebreak
    Where \c seedSampleSet is typically one of
    \plaintable
        uniform | uniform particle seeding
        cloud   | cloud of points
        patchSeed | seeding via patch faces
        triSurfaceMeshPointSet | points according to a tri-surface mesh
    \endplaintable

SeeAlso
    Foam::functionObject
    Foam::OutputFilterFunctionObject
    Foam::sampledSet
    Foam::streamLineBase
    Foam::streamLine

SourceFiles
    wallBoundedStreamLine.C

\*---------------------------------------------------------------------------*/

#ifndef wallBoundedStreamLine_H
#define wallBoundedStreamLine_H

#include "streamLineBase.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class wallBoundedStreamLine Declaration
\*---------------------------------------------------------------------------*/

class wallBoundedStreamLine
:
    public streamLineBase
{
   // Private Member Functions

        //- Find wall tet on cell
        tetIndices findNearestTet
        (
            const PackedBoolList& isWallPatch,
            const point& seedPt,
            const label cellI
        ) const;

        //- Disallow default bitwise copy construct
        wallBoundedStreamLine(const wallBoundedStreamLine&);

        //- Disallow default bitwise assignment
        void operator=(const wallBoundedStreamLine&);


public:

    //- Runtime type information
    TypeName("wallBoundedStreamLine");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        wallBoundedStreamLine
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~wallBoundedStreamLine();


    // Member Functions

        //- Read settings
        virtual void read(const dictionary&);

        //- Do the actual tracking to fill the track data
        virtual void track();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
