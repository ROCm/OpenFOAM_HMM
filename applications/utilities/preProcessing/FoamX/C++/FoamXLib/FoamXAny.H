/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    FoamXAny

Description
    Class to manage the mapping between Foam types, FoamX types and the
    CORBA Any type.

SourceFiles
    FoamXAny.C

\*---------------------------------------------------------------------------*/

#ifndef FoamX_H
#define FoamX_H

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

// Forward declaration of friend functions and operators

class FoamXAny;

void operator>
(
    const FoamXServer::FoamXAny& any1,
    const FoamXServer::FoamXAny& any2
);

void operator<
(
    const FoamXServer::FoamXAny& any1,
    const FoamXServer::FoamXAny& any2
);


/*---------------------------------------------------------------------------*\
                           Class Any Declaration
\*---------------------------------------------------------------------------*/

class FoamXAny
:
    public FoamXServer::FoamXAny
{
    // Private member functions

        void setDefaultValue();


public:

    // Constructors

        //- Default constructor
        FoamXAny(FoamXServer::FoamXType type = FoamXServer::Type_Undefined);

        //- Read constructor given type
        FoamXAny(FoamXServer::FoamXType type, Foam::Istream& is);

        //- Copy constructor
        FoamXAny(const FoamXAny& any);


    //- Destructor

        virtual ~FoamXAny();


    // Member Functions

        // Access

            void setType(FoamXServer::FoamXType type);

            //- Sets the value of this object. Throws an exception if the types
            //  are incompatible.
            void setValue(const FoamXServer::FoamXAny& any);

            //- Returns true if the value of this object has been set.
            bool IsSet() const
            {
                return (value.type()->kind() != CORBA::tk_null);
            }

            void setMin();
            void setMax();


        // Read

            void read(Foam::Istream& is);


        // Write

            void write(Foam::Ostream& os) const;


        // Member Operators

            void operator=(const FoamXServer::FoamXAny& any);
            void operator=(const FoamXAny& any);

            bool operator==(const FoamXServer::FoamXAny& any) const;
            bool operator!=(const FoamXServer::FoamXAny& any) const;

            friend void operator>
            (
                const FoamXServer::FoamXAny& any1,
                const FoamXServer::FoamXAny& any2
            );

            friend void operator<
            (
                const FoamXServer::FoamXAny& any1,
                const FoamXServer::FoamXAny& any2
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
