/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::functionObjects::surfaceDistance

Group
    grpFieldFunctionObjects

Description
    Writes a scalar field whose value is the distance to the nearest surface.
    The output field name is 'surfaceDistance'.

Usage
    Example of function object specification:
    \verbatim
    surfaceDistance
    {
        // Where to load it from
        libs            (libfieldFunctionObjects);
        type            surfaceDistance;

        geometry
        {
            motorBike.obj
            {
                type triSurfaceMesh;
                name motorBike;
            }
        }
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description                  | Required | Default
        type         | Type name: surfaceDistance   | yes |
        geometry     | Surfaces                     | yes  | no
        doCells      | Calculate distance from cell | no | true
    \endtable

    Run with e.g.

        postProcess -func surfaceDistance


See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::timeControl

SourceFiles
    surfaceDistance.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_surfaceDistance_H
#define functionObjects_surfaceDistance_H

#include "fvMeshFunctionObject.H"
#include "searchableSurfaces.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                  Class surfaceDistance Declaration
\*---------------------------------------------------------------------------*/

class surfaceDistance
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        //- Switch to calculate distance-to-cells
        Switch doCells_;

        //- Geometry
        autoPtr<searchableSurfaces> geomPtr_;


private:

    // Private member functions

        //- No copy construct
        surfaceDistance(const surfaceDistance&) = delete;

        //- No copy assignment
        void operator=(const surfaceDistance&) = delete;


public:

    //- Runtime type information
    TypeName("surfaceDistance");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        surfaceDistance
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~surfaceDistance();


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Calculate the interpolated fields
        virtual bool execute();

        //- Write the interpolated fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
