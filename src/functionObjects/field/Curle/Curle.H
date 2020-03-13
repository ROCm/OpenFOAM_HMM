/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::functionObjects::Curle

Group
    grpFieldFunctionObjects

Description
    Computes the acoustic pressure based on Curle's analogy.

    Curle's analogy is implemented as:

    \f[
        p' = \frac{1}{4 \pi c_0}\frac{\vec d}{|\vec d|^2}\cdot\frac{d\vec F}{dt}
    \f]

    where
    \vartable
      p'        | Curle's acoustic pressure     [Pa] or [Pa \f$(m^3/\rho)\f$]
      c_0       | Reference speed of sound              [m/s]
      \vec d    | Distance vector to observer locations [m]
      \vec F    | Force                         [N] or [N (\f$m^3/\rho\f$)]
    \endvartable

    Operands:
    \table
      Operand           | Type           | Location
      input             | volScalarField | $FOAM_CASE/\<time\>/\<inpField\>
      output file       | -              | -
      output field      | volScalarField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Note
    Only the normal-pressure force is included in the force calculation.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    Curle1
    {
        // Mandatory entries (unmodifiable)
        type            Curle;
        libs            (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        patches         (<patch1> <patch2> ... <patchN>)
        c0              343;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: Curle                   | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      patches      | Names of the operand patches       | wordList | yes | -
      c0           | Reference speed of sound [m/s]     | scalar | yes | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::fieldExpression
    - ExtendedCodeGuide::functionObjects::field::Curle

SourceFiles
    Curle.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_Curle_H
#define functionObjects_Curle_H

#include "fieldExpression.H"
#include "dimensionedScalar.H"
#include "dimensionedVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                            Class Curle Declaration
\*---------------------------------------------------------------------------*/

class Curle
:
    public fieldExpression
{
    // Private Data

        // Read from dictionary

            //- Patches to integrate forces over
            labelHashSet patchSet_;

            //- Area-averaged centre of patch faces
            dimensionedVector x0_;

            //- Reference speed of sound
            dimensionedScalar c0_;


protected:

<<<<<<< HEAD
    // Protected Member Functions

        //- Calculate acoustic pressure field and return true if successful
        virtual bool calc();

        //- No copy construct
        Curle(const Curle&) = delete;

        //- No copy assignment
        void operator=(const Curle&) = delete;

=======
        //- Calculate the acoustic pressure field and return true if successful
        virtual bool calc();

>>>>>>> DOC: elaborate the usage of function objects

public:

    //- Runtime type information
    TypeName("Curle");


    // Constructors

        //- Construct from Time and dictionary
        Curle
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        Curle(const Curle&) = delete;

        //- No copy assignment
        void operator=(const Curle&) = delete;


    //- Destructor
    virtual ~Curle() = default;


    // Member Functions

        //- Read the Curle data
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
