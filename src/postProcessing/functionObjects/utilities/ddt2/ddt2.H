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
    Foam::ddt2

Group
    grpFVFunctionObjects

Description
    This function object calculates the magnitude squared
    of d(scalarField)/dt.

    The result can be used further for determining variance or RMS values
    (for example).

    Example of function object specification:
    \verbatim
    dpdt2
    {
        type        ddt2;
        functionObjectLibs ("libutilityFunctionObjects.so");
        fields      (p);
        result      d@@dt2;
        ...
    }
    \endverbatim

    \heading Function object usage
    \table
        Property | Description                | Required  | Default value
        type     | type name: ddt2            | yes       |
        fields   | Name of fields to process  | yes       |
        result   | Name of results            | no        | magSqr(ddt(@@))
        log      | Log to standard output     | no        | yes
        mag      | Use 'mag' instead of 'magSqr' | no     | false
    \endtable

    Note that the optional 'mag' entry cannot be changed during the simulation
    since it alters the dimensions of the output field.

    A list of fields can contain exact names or regular expressions.
    The token '\@\@' in the result name is replaced by the name of the source
    field.

    The function object will skip over fields that appear to have
    already been processed (ie, their names are similar to the output names).

SourceFiles
    ddt2.C
    IOddt2.H

\*---------------------------------------------------------------------------*/

#ifndef ddt2_H
#define ddt2_H

#include "volFieldsFwd.H"
#include "OFstream.H"
#include "wordReList.H"
#include "regExp.H"
#include "HashSet.H"
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
                            Class ddt2 Declaration
\*---------------------------------------------------------------------------*/

class ddt2
{
    // Private data

        //- Name of this ddt2 object
        word name_;

        //- Reference to the database
        const objectRegistry& obr_;

        //- On/off switch
        bool active_;

        //- Name of fields to process
        wordReList selectFields_;

        //- Formatting for the result fields.
        word resultName_;

        //- Avoid processing the same field twice.
        mutable regExp blacklist_;

        //- Hashed names of result fields.
        wordHashSet results_;

        //- Switch to send output to Info as well as to file
        Switch log_;

        //- Use 'mag' instead of 'magSqr'.
        //  Cannot be adjusted during the simulation since it alters the
        //  dimensions of the output field.
        const Switch mag_;


    // Private Member Functions

        //- Check that the word contains the appropriate substitution token(s).
        static bool checkFormatName(const word&);

        //- Eliminate duplicate 'word' entries
        static void uniqWords(wordReList&);


        //- Accept unless field name appears to have already been processed
        bool accept(const word& fieldName) const;

        //- Apply for the volume field type
        template<class FieldType>
        int apply(const fvMesh& mesh, const word& inputName, int& state);

        //- Process by trying to apply for various volume field types.
        int process(const fvMesh& mesh, const word& inputName);

        //- Calculate the ddt2 fields
        void process(const fvMesh& mesh);


        //- Disallow default bitwise copy construct
        ddt2(const ddt2&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const ddt2&) = delete;


public:

    //- Runtime type information
    TypeName("ddt2");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        ddt2
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~ddt2();


    // Member Functions

        //- Return name of the ddt2 function object
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the ddt2 specification
        virtual void read(const dictionary&);

        //- Calculate the ddt2 fields
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end()
        {}

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet()
        {}

        //- Write the ddt fields
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
