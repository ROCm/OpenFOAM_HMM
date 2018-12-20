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
    Foam::functionObjects::fieldExtents

Group
    grpFieldFunctionObjects

Description
    Calculates the spatial minimum and maximum extents of a field

    The extents are derived from the bound box limits after identifying the
    locations where field values exceed the user-supplied threshold value.

Usage
    Example of function object specification:
    \verbatim
    fieldExtents1
    {
        type        fieldExtents;
        libs        ("libfieldFunctionObjects.so");
        ...
        writeToFile yes;
        log         yes;
        fields      (alpha);
        threshold   0.5;
        patches     ();
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property      | Description              | Required   | Default value
        type          | type name: fieldExtents  | yes        |
        writeToFile   | write extents data to file | no       | yes
        log           | write extents data to standard output | no | yes
        internalField | Process the internal field | no       | yes
        threshold     | Field value to identify extents boundary | yes |
        referencePosition | Reference position   | no         | (0 0 0)
        fields        | list of fields to process | yes       |
        patches       | list of patches to process | no       | \<all patches\>
    \endtable

    Output data is written to the file \<timeDir\>/fieldExtents.dat

Note
    For non-scalar fields, the magnitude of the field is employed and compared
    to the threshold value.

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    fieldExtents.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fieldExtents_H
#define functionObjects_fieldExtents_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class fieldExtents Declaration
\*---------------------------------------------------------------------------*/

class fieldExtents
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data

        //- Flag to write the internal field extents
        bool internalField_;

        //- Threshold value
        scalar threshold_;

        //- Reference position; default = (0 0 0)
        point C0_;

        //- Fields to assess
        volFieldSelection fieldSet_;

        //- Patches to assess
        labelHashSet patchIDs_;


    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);

        //- Return the field mask
        template<class Type>
        tmp<volScalarField> calcMask
        (
            const GeometricField<Type, fvPatchField, volMesh>& field
        ) const;

        //- Main calculation
        template<class Type>
        void calcFieldExtents
        (
            const word& fieldName,
            const bool calcMag = false
        );

        //- No copy construct
        fieldExtents(const fieldExtents&) = delete;

        //- No copy assignment
        void operator=(const fieldExtents&) = delete;


public:

    //- Runtime type information
    TypeName("fieldExtents");


    // Constructors

        //- Construct from Time and dictionary
        fieldExtents
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fieldExtents() = default;


    // Member Functions

        //- Read the field extents data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the fieldExtents
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<volScalarField> fieldExtents::calcMask
(
    const GeometricField<scalar, fvPatchField, volMesh>& field
) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldExtentsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
