/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::basicSolidThermo

Description
    The thermophysical properties of a basicSolidThermo

SourceFiles
    basicSolidThermo.C
    newBasicSolidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef basicSolidThermo_H
#define basicSolidThermo_H

#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class basicSolidThermo Declaration
\*---------------------------------------------------------------------------*/

class basicSolidThermo
:
    public IOdictionary
{
protected:

    // Protected data

        const fvMesh& mesh_;

        //- Temperature [K]
        volScalarField T_;


public:

    //- Runtime type information
    TypeName("basicSolidThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicSolidThermo,
        mesh,
        (const fvMesh& mesh),
        (mesh)
    );


    // Constructors

        //- Construct from mesh
        basicSolidThermo(const fvMesh&);

        //- Return a pointer to a new basicSolidThermo created from
        //  the solidThermophysicalProperties dictionary
        static autoPtr<basicSolidThermo> New(const fvMesh&);


    //- Destructor
    virtual ~basicSolidThermo();


    // Member Functions

        //- Update properties
        virtual void correct() = 0;


        // Physical constants which define the basicSolidThermo

            //- Temperature [K]
            inline const volScalarField& T() const;

            //- Temperature [K]
            inline volScalarField& T();

            //- Density [kg/m3]
            virtual tmp<volScalarField> rho() const = 0;

            //- Specific heat capacity [J/(kg.K)]
            virtual tmp<volScalarField> cp() const = 0;

            //- Thermal conductivity [W/(m.K)]
            virtual tmp<volScalarField> K() const = 0;

            //- Thermal conductivity [W/(m.K)]
            virtual tmp<volSymmTensorField> directionalK() const = 0;

            //- Heat of formation [J/kg]
            virtual tmp<volScalarField> Hf() const = 0;

            //- Emissivity []
            virtual tmp<volScalarField> emissivity() const = 0;


            // Per patch calculation

                //- Density [kg/m3]
                virtual tmp<scalarField> rho(const label patchI) const = 0;

                //- Specific heat capacity [J/(kg.K)]
                virtual tmp<scalarField> cp(const label patchI) const = 0;

                //- Thermal conductivity [W/(m.K)]
                virtual tmp<scalarField> K(const label patchI) const = 0;

                //- Thermal conductivity [W/(m.K)]
                virtual tmp<symmTensorField> directionalK(const label) const =0;

                //- Heat of formation [J/kg]
                virtual tmp<scalarField> Hf(const label patchI) const = 0;

                //- Emissivity []
                virtual tmp<scalarField> emissivity(const label) const = 0;


//        // Point wise properties
//
//            //- Density [kg/m3]
//            virtual scalar rho(const scalar T) const = 0;
//
//            //- Specific heat capacity [J/(kg.K)]
//            virtual scalar cp(const scalar T) const = 0;
//
//            //- Thermal conductivity [W/(m.K)]
//            virtual scalar K(const scalar T) const = 0;
//
//            //- Heat of formation [J/kg]
//            virtual scalar Hf(const scalar T) const = 0;
//
//            //- Emissivity []
//            virtual scalar emissivity(const scalar T) const = 0;


        // I-O

            //- Write the basicSolidThermo properties
            virtual bool writeData(Ostream& os) const = 0;

            //- Read solidThermophysicalProperties dictionary
            virtual bool read() = 0;


        // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const basicSolidThermo& s);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicSolidThermoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
