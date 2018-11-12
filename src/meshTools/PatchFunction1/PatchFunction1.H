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
    Foam::PatchFunction1

Description
    Top level data entry class for use in dictionaries. Provides a mechanism
    to specify a variable as a certain type, e.g. constant or time varying, and
    provide functions to return the (interpolated) value, and integral between
    limits.

    Extends the Function1 class by adding autoMap and rMap functions

SourceFiles
    PatchFunction1.C
    PatchFunction1New.C

SeeAlso
    Foam::Function1

\*---------------------------------------------------------------------------*/

#ifndef PatchFunction1_H
#define PatchFunction1_H

#include "dictionary.H"
#include "Field.H"
#include "polyPatch.H"
#include "coordinateScaling.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class Time;

// Forward declaration of friend functions and operators
template<class Type>
class PatchFunction1;

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const PatchFunction1<Type>&
);

/*---------------------------------------------------------------------------*\
                       Class PatchFunction1 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class PatchFunction1
:
    public refCount
{
    // Private Member Functions

        //- No copy assignment
        void operator=(const PatchFunction1<Type>&) = delete;


protected:

    // Protected data

        //- Name of entry
        const word name_;

        //- Reference to the patch
        const polyPatch& patch_;

        //- Whether to generate face or point values on patch
        const bool faceValues_;

        //- Optional local co-ordinate system and scaling
        coordinateScaling<Type> coordSys_;


public:

    typedef Field<Type> returnType;

    //- Runtime type information
    TypeName("PatchFunction1")

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        PatchFunction1,
        dictionary,
        (
            const polyPatch& pp,
            const word& type,
            const word& entryName,
            const dictionary& dict,
            const bool faceValues
        ),
        (pp, type, entryName, dict, faceValues)
    );


    // Constructors

        //- Construct from polyPatch and entry name
        PatchFunction1
        (
            const polyPatch& pp,
            const word& entryName,
            const bool faceValues = true
        );

        //- Construct from polyPatch, dictionary and entry name
        PatchFunction1
        (
            const polyPatch& pp,
            const word& entryName,
            const dictionary& dict,
            const bool faceValues = true
        );

        //- Copy constructor
        explicit PatchFunction1(const PatchFunction1<Type>& pf1);

        //- Copy constructor setting patch
        explicit PatchFunction1
        (
            const PatchFunction1<Type>& pf1,
            const polyPatch& pp
        );

        //- Construct and return a clone
        virtual tmp<PatchFunction1<Type>> clone() const = 0;

        //- Construct and return a clone setting patch
        virtual tmp<PatchFunction1<Type>> clone(const polyPatch& pp) const = 0;


    //- Selector
    static autoPtr<PatchFunction1<Type>> New
    (
        const polyPatch& pp,
        const word& entryName,
        const dictionary& dict,
        const bool faceValues = true
    );


    //- Destructor
    virtual ~PatchFunction1() = default;


    // Member Functions

        // Access

            //- Return the name of the entry
            const word& name() const;


        // Manipulation

            //- Convert time
            virtual void convertTimeBase(const Time& t);


        // Evaluation

            //- Return value as a function of (scalar) independent variable
            virtual tmp<Field<Type>> value(const scalar x) const;

            //- Is value constant (i.e. independent of x)
            virtual bool constant() const = 0;

            //- Is value uniform (i.e. independent of coordinate)
            virtual bool uniform() const = 0;

            //- Integrate between two (scalar) values
            virtual tmp<Field<Type>> integrate
            (
                const scalar x1,
                const scalar x2
            ) const;

            //- Helper: optionally convert coordinates to local coordinates
            virtual tmp<pointField> localPosition
            (
                const pointField& globalPos
            ) const;

            //- Apply optional transformation
            virtual tmp<Field<Type>> transform(const Field<Type>& fld) const;

            //- Apply optional transformation
            virtual tmp<Field<Type>> transform
            (
                const tmp<Field<Type>>& tfld
            ) const;


        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const FieldMapper& mapper);

            //- Reverse map the given PatchFunction1 onto this PatchFunction1
            virtual void rmap
            (
                const PatchFunction1<Type>& pf1,
                const labelList& addr
            );


        // I/O

            //- Ostream Operator
            friend Ostream& operator<< <Type>
            (
                Ostream& os,
                const PatchFunction1<Type>& func
            );

            //- Write in dictionary format
            virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makePatchFunction1(Type)                                               \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(PatchFunction1<Type>, 0);              \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        PatchFunction1<Type>,                                                  \
        dictionary                                                             \
    );

#define makePatchFunction1Type(SS, Type)                                       \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(PatchFunction1Types::SS<Type>, 0);     \
                                                                               \
    PatchFunction1<Type>::adddictionaryConstructorToTable                      \
        <PatchFunction1Types::SS<Type>>                                        \
        add##SS##Type##ConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PatchFunction1.C"
    #include "PatchFunction1New.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
