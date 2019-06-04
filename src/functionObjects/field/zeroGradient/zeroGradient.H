/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::functionObjects::zeroGradient

Group
    grpFieldFunctionObjects

Description
    This function object creates a volume field with zero-gradient
    boundary conditions from another volume field.

    The result can be used, for example, to post-process near-wall
    field values.

Usage
    Example of function object specification:
    \verbatim
    zeroGrad
    {
        type        zeroGradient;
        libs        ("libfieldFunctionObjects.so");
        fields      (U "(T|k|epsilon|omega)");
        result      @@nearWall;
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property | Description                | Required  | Default value
        type     | type name: zeroGradient    | yes       |
        fields   | Name of fields to process  | yes       |
        result   | Name of results            | no        | zeroGradient(@@)
        log      | Log to standard output     | no        | no
    \endtable

    A list of fields can contain exact names or regular expressions.
    The token '\@\@' in the result name is replaced by the name of the source
    field. In the special case of a single source field (specified as
    a non-regex), the '\@\@' token checking is suppressed.

    The function object will skip over fields that would not benefit
    - ie, only processor, empty, zeroGradient, symmetry patches.
    This check should also prevent processing fields multiple times.

SourceFiles
    zeroGradient.C
    zeroGradientTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_zeroGradient_H
#define functionObjects_zeroGradient_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"
#include "OFstream.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class zeroGradient Declaration
\*---------------------------------------------------------------------------*/

class zeroGradient
:
    public fvMeshFunctionObject
{
    // Private data

        //- Name of fields to process.
        wordRes selectFields_;

        //- Formatting for the result fields.
        word resultName_;

        //- Hashed names of result fields, and their type.
        HashTable<word> results_;


    // Private Member Functions

        //- Check that string contains the appropriate substitution token(s).
        static bool checkFormatName(const std::string& str);


        //- Accept unless field only has constraint patches
        //  (ie, empty/zero-gradient/processor).
        //  This should also avoid fields that were already processed by
        //  zeroGradient.
        template<class Type>
        static bool accept(const GeometricField<Type, fvPatchField, volMesh>&);

        //- Apply for the volume field type
        template<class Type>
        int apply(const word& inputName, int& state);

        //- Process by trying to apply for various volume field types.
        int process(const word& inputName);


        //- No copy construct
        zeroGradient(const zeroGradient&) = delete;

        //- No copy assignment
        void operator=(const zeroGradient&) = delete;


public:

    //- Runtime type information
    TypeName("zeroGradient");


    // Constructors

        //- Construct from Time and dictionary
        zeroGradient
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~zeroGradient() = default;


    // Member Functions

        //- Read the zeroGradient specification
        virtual bool read(const dictionary& dict);

        //- Calculate the zeroGradient fields
        virtual bool execute();

        //- Write the zeroGradient fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "zeroGradientTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
