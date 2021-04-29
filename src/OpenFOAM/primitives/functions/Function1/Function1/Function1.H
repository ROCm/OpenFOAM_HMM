/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::Function1

Description
    Top level data entry class for use in dictionaries. Provides a mechanism
    to specify a variable as a certain type, e.g. constant or table, and
    provide functions to return the (interpolated) value, and integral between
    limits.

    The New factory method attempts to deal with varying types of input.
    It accepts primitive or dictionary entries for dispatching to different
    function types, but wraps unspecified types as "constant".

    In the dictionary form, the coefficents are the dictionary itself.
    This is arguably the more readable form.
    For example,
    \verbatim
        <entryName>
        {
            type      linearRamp;
            start     10;
            duration  20;
        }
    \endverbatim

    In the primitive form, the coefficents are provided separately.
    For example,
    \verbatim
        <entryName>   linearRamp;
        <entryName>Coeffs
        {
            start     10;
            duration  20;
        }
    \endverbatim
    The coeffs dictionary is optional, since it is not required by all types.
    For example,
    \verbatim
        <entryName>   zero;
    \endverbatim

SourceFiles
    Function1.C
    Function1New.C

\*---------------------------------------------------------------------------*/

#ifndef Function1_H
#define Function1_H

#include "function1Base.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;
template<class Type> class Function1;
template<class Type> Ostream& operator<<(Ostream&, const Function1<Type>&);

/*---------------------------------------------------------------------------*\
                          Class Function1 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Function1
:
    public function1Base
{
    // Private Member Functions

        //- Selector, with alternative entry, fallback redirection, etc
        static autoPtr<Function1<Type>> New
        (
            const word& entryName,   // Entry name for function
            const entry* eptr,       // Eg, dict.findEntry(entryName)
            const dictionary& dict,
            const word& redirectType,
            const bool mandatory
        );


protected:

    // Protected Member Functions

        //- No copy assignment
        void operator=(const Function1<Type>&) = delete;


public:

    typedef Type returnType;

    //- Runtime type information
    TypeName("Function1")

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        Function1,
        dictionary,
        (
            const word& entryName,
            const dictionary& dict
        ),
        (entryName, dict)
    );


    // Constructors

        //- Construct from entry name
        explicit Function1(const word& entryName);

        //- Construct from entry name and dictionary (unused)
        Function1(const word& entryName, const dictionary& dict);

        //- Copy construct
        explicit Function1(const Function1<Type>& rhs);

        //- Construct and return a clone
        virtual tmp<Function1<Type>> clone() const = 0;


    // Selectors

        //- Selector, with fallback redirection
        static autoPtr<Function1<Type>> New
        (
            const word& entryName,
            const dictionary& dict,
            const word& redirectType,
            const bool mandatory = true
        );

        //- Compatibility selector, with fallback redirection
        static autoPtr<Function1<Type>> NewCompat
        (
            const word& entryName,
            std::initializer_list<std::pair<const char*,int>> compat,
            const dictionary& dict,
            const word& redirectType = word::null,
            const bool mandatory = true
        );

        //- Selector, without fallback redirection
        static autoPtr<Function1<Type>> New
        (
            const word& entryName,
            const dictionary& dict,
            const bool mandatory = true
        );

        //- An optional selector
        static autoPtr<Function1<Type>> NewIfPresent
        (
            const word& entryName,
            const dictionary& dict,
            const word& redirectType = word::null
        );


    //- Destructor
    virtual ~Function1() = default;


    // Member Functions

        // Evaluation

            //- Return value as a function of (scalar) independent variable
            virtual Type value(const scalar x) const;

            //- Return value as a function of (scalar) independent variable
            virtual tmp<Field<Type>> value(const scalarField& x) const;

            //- Integrate between two (scalar) values
            virtual Type integrate(const scalar x1, const scalar x2) const;

            //- Integrate between two (scalar) values
            virtual tmp<Field<Type>> integrate
            (
                const scalarField& x1,
                const scalarField& x2
            ) const;


    // I/O

        //- Ostream Operator
        friend Ostream& operator<< <Type>
        (
            Ostream& os,
            const Function1<Type>& func
        );


        //- Write in dictionary format.
        //  \note The base output is \em without an END_STATEMENT
        virtual void writeData(Ostream& os) const;

        //- Write coefficient entries in dictionary format
        void writeEntries(Ostream& os) const;
};


/*---------------------------------------------------------------------------*\
                        Class FieldFunction1 Declaration
\*---------------------------------------------------------------------------*/

template<class Function1Type>
class FieldFunction1
:
    public Function1Type
{
public:

    typedef typename Function1Type::returnType Type;


    // Constructors

        //- Construct from entry name and dictionary
        FieldFunction1(const word& entryName, const dictionary& dict);

        //- Construct and return a clone
        virtual tmp<Function1<Type>> clone() const;


    //- Destructor
    virtual ~FieldFunction1() = default;


    // Member Functions

        // Evaluation

            using Function1Type::value;
            using Function1Type::integrate;

            //- Return value as a function of (scalar) independent variable
            virtual tmp<Field<Type>> value(const scalarField& x) const;

            //- Integrate between two (scalar) values
            virtual tmp<Field<Type>> integrate
            (
                const scalarField& x1,
                const scalarField& x2
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Define Function1 run-time selection
#define makeFunction1(Type)                                                    \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(Function1<Type>, 0);                   \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        Function1<Type>,                                                       \
        dictionary                                                             \
    );


// Define (templated) Function1, add to (templated) run-time selection
#define makeFunction1Type(SS, Type)                                            \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(Function1Types::SS<Type>, 0);          \
                                                                               \
    Function1<Type>::adddictionaryConstructorToTable                           \
        <FieldFunction1<Function1Types::SS<Type>>>                             \
        add##SS##Type##ConstructorToTable_;


// Define a non-templated Function1 and add to (templated) run-time selection
#define makeConcreteFunction1(SS, Type)                                        \
                                                                               \
    defineTypeNameAndDebug(SS, 0);                                             \
                                                                               \
    Function1<Type>::adddictionaryConstructorToTable                           \
        <FieldFunction1<SS>>                                                   \
        add##SS##Type##ConstructorToTable_;


// Define scalar Function1 and add to (templated) run-time selection
#define makeScalarFunction1(SS)                                                \
                                                                               \
    makeConcreteFunction1(SS, scalar);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Function1.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
