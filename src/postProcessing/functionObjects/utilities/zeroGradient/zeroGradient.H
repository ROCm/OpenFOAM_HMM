/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::zeroGradient

Group
    grpFVFunctionObjects

Description
    This function object creates a volume field with zero-gradient
    boundary conditions from another volume field.

    The result can be used, for example, to post-process near-wall
    field values.

    Example of function object specification:
    \verbatim
    zeroGrad
    {
        type        zeroGradient;
        functionObjectLibs ("libutilityFunctionObjects.so");
        fields      (U "(T|k|epsilon|omega)");
        result      @@nearWall;
        ...
    }
    \endverbatim

    \heading Function object usage
    \table
        Property | Description                | Required  | Default value
        type     | type name: zeroGradient    | yes       |
        fields   | Name of fields to process  | yes       |
        result   | Name of results            | no        | zeroGradient(@@)
        log      | Log to standard output     | no        | no
    \endtable

    A list of fields can contain exact names or regular expressions.
    The token '\@\@' in the result name is replaced by the name of the source
    field.

    The function object will skip over fields that would not benefit
    - ie, only processor, empty, zeroGradient, symmetry patches.
    This check should also prevent processing fields multiple times.

SourceFiles
    zeroGradient.C
    zeroGradientFunctionObject.C
    IOzeroGradient.H

\*---------------------------------------------------------------------------*/

#ifndef zeroGradient_H
#define zeroGradient_H

#include "volFieldsFwd.H"
#include "OFstream.H"
#include "wordReList.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class fvMesh;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                            Class zeroGradient Declaration
\*---------------------------------------------------------------------------*/

class zeroGradient
{
    // Private data

        //- Name of this zeroGradient object
        word name_;

        //- Reference to the database
        const objectRegistry& obr_;

        //- On/off switch
        bool active_;

        //- Name of fields to process
        wordReList selectFields_;

        //- Formatting for the result fields.
        word resultName_;

        //- Hashed names of result fields, and their type
        HashTable<word> results_;

        //- Switch to send output to Info as well as to file
        Switch log_;


    // Private Member Functions

        //- Check that the word contains the appropriate substitution token(s).
        static bool checkFormatName(const word&);

        //- Eliminate duplicate 'word' entries
        static void uniqWords(wordReList&);


        //- Accept unless field only has empty/zero-gradient/processor patches
        template<class Type>
        static bool accept(const GeometricField<Type, fvPatchField, volMesh>&);

        //- Apply for the volume field type
        template<class Type>
        int apply(const fvMesh& mesh, const word& inputName, int& state);

        //- Process by trying to apply for various volume field types.
        int process(const fvMesh& mesh, const word& inputName);

        //- Calculate the zeroGradient fields
        void process(const fvMesh& mesh);


        //- Disallow default bitwise copy construct
        zeroGradient(const zeroGradient&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const zeroGradient&) = delete;


public:

    //- Runtime type information
    TypeName("zeroGradient");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        zeroGradient
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~zeroGradient();


    // Member Functions

        //- Return name of the zeroGradient function object
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the zeroGradient specification
        virtual void read(const dictionary&);

        //- Calculate the zeroGradient fields
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end()
        {}

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet()
        {}

        //- Write the zeroGradient fields
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
