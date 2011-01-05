/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2011 OpenCFD Ltd.
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
    Foam::pointSolid

Description
    The thermophysical properties of a pointSolid

SourceFiles
    pointSolid.C

\*---------------------------------------------------------------------------*/

#ifndef pointSolid_H
#define pointSolid_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointSolid;

Ostream& operator<<
(
    Ostream&,
    const pointSolid&
);


/*---------------------------------------------------------------------------*\
                       Class pointSolid Declaration
\*---------------------------------------------------------------------------*/

class pointSolid
{

    // Private data

        //- Density [kg/m3]
        scalar rho_;

        //- Specific heat capacity [J/(kg.K)]
        scalar Cp_;

        //- Thermal conductivity [W/(m.K)]
        scalar K_;

        //- Heat of formation [J/kg]
        scalar Hf_;

        //- Emissivity
        scalar emissivity_;


public:

    //- Runtime type information
    TypeName("pointSolid");


    // Declare run-time constructor selection tables

    declareRunTimeSelectionTable
    (
        autoPtr,
        pointSolid,
        ,
        (),
        ()
    );

    declareRunTimeSelectionTable
    (
        autoPtr,
        pointSolid,
        Istream,
        (Istream& is),
        (is)
    );

    declareRunTimeSelectionTable
    (
        autoPtr,
        pointSolid,
        dictionary,
        (const dictionary& dict),
        (dict)
    );


    // Constructors

        //- Construct from components
        pointSolid
        (
            scalar rho,
            scalar Cp,
            scalar K,
            scalar Hf,
            scalar emissivity
        );

        //- Construct from Istream
        pointSolid(Istream& is);

        //- Construct from dictionary
        pointSolid(const dictionary& dict);

        //- Construct copy
        pointSolid(const pointSolid& s);

        //- Construct and return clone
        virtual autoPtr<pointSolid> clone() const
        {
            return autoPtr<pointSolid>(new pointSolid(*this));
        }


    // Selectors

        //- Return a pointer to a new pointSolid created from input
        static autoPtr<pointSolid> New(Istream& is);

        //- Return a pointer to a new pointSolid created from dictionary
        static autoPtr<pointSolid> New(const dictionary& dict);


    //- Destructor
    virtual ~pointSolid()
    {}


    // Member Functions

        // Physical constants which define the pointSolid

            //- Density [kg/m3]
            inline scalar rho() const;

            //- Specific heat capacity [J/(kg.K)]
            inline scalar Cp() const;

            //- Thermal conductivity [W/(m.K)]
            inline scalar K() const;

            //- Heat of formation [J/kg]
            inline scalar Hf() const;

            //- Total enthalpy - reference to Tstd [J/kg]
            inline scalar H(const scalar T) const;

            //- Sensible enthalpy - reference to Tstd [J/kg]
            inline scalar Hs(const scalar T) const;

            //- Emissivity []
            inline scalar emissivity() const;


        // I-O

            //- Write the pointSolid properties
            virtual void writeData(Ostream& os) const;


        // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const pointSolid& s);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pointSolidI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
