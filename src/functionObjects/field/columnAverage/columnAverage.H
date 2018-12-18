/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

Class
    Foam::functionObjects::columnAverage

Group
    grpFieldFunctionObjects

Description
    Averages columns of cells for layered meshes.

    For each patch face, calculates the average value of all cells attached in
    the patch face normal direction, and then pushes the average value back
    to all cells in the column.

    Useful for channel-like cases where we want to average fields in the
    spanwise direction.

Usage
    Example of function object specification:
    \verbatim
    columnAverage1
    {
        type        columnAverage;
        libs        ("libfieldFunctionObjects.so");
        ...
        patches     (front side);
        fields      (U p);
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description               | Required    | Default value
        type         | type name: fieldMinMax    | yes         |
        patches      | list of patches to collapse onto | yes  |
        fields       | list of fields to process | yes         |
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    columnAverage.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_columnAverage_H
#define functionObjects_columnAverage_H

#include "volFieldsFwd.H"
#include "fvMeshFunctionObject.H"
#include "volFieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class globalIndex;
class meshStructure;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class columnAverage Declaration
\*---------------------------------------------------------------------------*/

class columnAverage
:
    public fvMeshFunctionObject
{
    // Private data

        //- Patches on which to collapse the fields
        labelHashSet patchSet_;

        //- Fields to collapse
        volFieldSelection fieldSet_;

        mutable autoPtr<globalIndex> globalFaces_;
        mutable autoPtr<globalIndex> globalEdges_;
        mutable autoPtr<globalIndex> globalPoints_;
        mutable autoPtr<meshStructure> meshStructurePtr_;


    // Private Member Functions

        //- Create the column average field name
        const word averageName(const word& fieldName) const;

        //- Return the column-based addressing
        const meshStructure& meshAddressing(const polyMesh&) const;

        //- Calculate the averaged field and return true if successful
        template<class Type>
        bool columnAverageField(const word& fieldName);


public:

    //- Runtime type information
    TypeName("columnAverage");


    // Constructors

        //- Construct from Time and dictionary
        columnAverage
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~columnAverage() = default;


    // Member Functions

        //- Read the settings
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the results
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "columnAverageTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
