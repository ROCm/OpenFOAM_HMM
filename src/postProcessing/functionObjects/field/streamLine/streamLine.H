/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::streamLine

Group
    grpFieldFunctionObjects

Description
    This function object generates streamline data by sampling a set of
    user-specified fields along a particle track, transported by a
    user-specified velocity field.

    Example of function object specification:
    \verbatim
    streamLine1
    {
        type            streamLine;
        functionObjectLibs ("libfieldFunctionObjects.so");
        ...
        setFormat       vtk;
        UName           U;
        trackForward    yes;
        fields
        (
            U
            p
        );
        lifeTime        10000;
        trackLength     1e-3;
        nSubCycle       5;
        bounds          (0.2 -10 -10)(0.22 10 10);
        cloudName       particleTracks;
        seedSampleSet   uniform;
        uniformCoeffs
        {
            type        uniform;
            axis        x;  //distance;
            start       (-0.0205 0.0001 0.00001);
            end         (-0.0205 0.0005 0.00001);
            nPoints     100;
        }
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: streamLine   | yes         |
        setFormat    | output data type        | yes         |
        UName        | tracking velocity field name | yes    |
        fields       | fields to sample        | yes         |
        lifetime     | maximum number of particle tracking steps | yes |
        trackLength  | tracking segment length | no          |
        nSubCycle    | number of tracking steps per cell | no|
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
        triSurfaceMeshPointSet | points according to a tri-surface mesh
    \endplaintable

Note
    When specifying the track resolution, the \c trackLength OR \c nSubCycle
    option should be used

SeeAlso
    Foam::functionObject
    Foam::OutputFilterFunctionObject
    Foam::sampledSet
    Foam::wallBoundedStreamLine
    Foam::streamLineBase

SourceFiles
    streamLine.C

\*---------------------------------------------------------------------------*/

#ifndef streamLine_H
#define streamLine_H

#include "streamLineBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;

/*---------------------------------------------------------------------------*\
                         Class streamLine Declaration
\*---------------------------------------------------------------------------*/

class streamLine
:
    public streamLineBase
{
    // Private data

        //- Number of subcycling steps
        label nSubCycle_;


   // Private Member Functions

        //- Disallow default bitwise copy construct
        streamLine(const streamLine&);

        //- Disallow default bitwise assignment
        void operator=(const streamLine&);


public:

    //- Runtime type information
    TypeName("streamLine");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        streamLine
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~streamLine();


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
