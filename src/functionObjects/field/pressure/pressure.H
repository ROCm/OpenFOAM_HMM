/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::functionObjects::pressure

Group
    grpFieldFunctionObjects

Description
    Includes tools to manipulate the pressure into different forms.

    These currently include:
    - static pressure
        \f[
            p_s = p_{ref} + \rho p_k
        \f]
    - total pressure
        \f[
            p_0 = p_{ref} + p + 0.5 \rho |U|^2
        \f]
    - isentropic pressure
        \f[
            p_i = p*(1 + ((gamma-1)*M^2)/2)^(gamma/(gamma - 1))
        \f]
    - static pressure coefficient
        \f[
            Cp = \frac{p_s - p_{\inf}}{0.5 \rho_{\inf} |U_{\inf}|^2}
        \f]
    - total pressure coefficient
        \f[
            Cp_0 = \frac{p_0 - p_{\inf}}{0.5 \rho_{\inf} |U_{\inf}|^2}
        \f]

    where
    \vartable
        \rho        | Density [kg/m3]
        U           | Velocity [m/s]
        \rho_{\inf} | Freestream density [kg/m3]
        p_{\inf}    | Freestream pressure [Pa]
        U_{\inf}    | Freestream velocity [m/s]
        p_k         | Kinematic pressure (p/rho)[m2/s2]
        p_s         | Statoc pressure [Pa]
        p_0         | Total pressure [Pa]
        p_{ref}     | Reference pressure level [Pa]
        p_i         | Total isentropic pressure
        Cp          | Pressure coefficient
        Cp_0        | Total pressure coefficient
    \endvartable

    The function object will operate on both kinematic (\f$ p_k \f$) and static
    pressure (\f$ p \f$) fields, and the result is written as a
    volScalarField.

Usage
    Example of function object specification to calculate pressure coefficient:
    \verbatim
    pressure1
    {
        type        pressure;
        libs        ("libfieldFunctionObjects.so");
        ...
        mode        staticCoeff;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description                 | Required   | Default value
        type         | type name: pressure         | yes        |
        field        | Name of the pressure field  | no         | p
        U            | Name of the velocity field  | no         | U
        rho          | Name of the density field   | no         | rho
        result       | Name of the resulting field | no         | derived from p
        mode         | Calculation mode (see below) | yes       |
        pRef         | Reference pressure for total pressure | no | 0
        pInf         | Freestream pressure for coefficient calculation | no |
        UInf         | Freestream velocity for coefficient calculation | no |
        rhoInf       | Freestream density for coefficient calculation  | no |
        hydrostaticMode | Hydrostatic contributions (see below) | no | none
        g            | Gravity vector (see below) | no |
        hRef         | Reference height (see below) | no |
    \endtable

    The \c mode entry is used to select the type of pressure that is calculated.
    Selections include:
    - static
    - total
    - isentropic
    - staticCoeff
    - totalCoeff

    The optional \c hydrostaticMode entry provides handling for the term
    \f$ \rho (\vec{g} \dot \vec{h})\f$ where options include
    - \c none : not included
    - \c add : add the term, e.g. to convert from p_rgh to p
    - \c subtract : subtract the term, e.g. to convert from p to p_rgh

    If the \c hydrostaticMode is active, values are also required for
    gravity, \c g, and reference height, \c hRef.  By default these will be
    retrieved from the database.  When these values are not available
    the user must provide them, e.g.
    \verbatim
        g           (0 -9.81 0);
        hRef        0;
    \endverbatim


See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    pressure.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_pressure_H
#define functionObjects_pressure_H

#include "fieldExpression.H"
#include "volFieldsFwd.H"
#include "dimensionedVector.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class pressure Declaration
\*---------------------------------------------------------------------------*/

class pressure
:
    public fieldExpression
{
public:

    // Public Data Types

        //- Enumeration for pressure calculation mode
        enum mode : unsigned
        {
            STATIC       = (1 << 0), //!< Static pressure
            TOTAL        = (1 << 1), //!< Total pressure
            ISENTROPIC   = (1 << 2), //!< Isentropic pressure
            COEFF        = (1 << 3), //!< Coefficient manipulator
            STATIC_COEFF = (STATIC | COEFF),
            TOTAL_COEFF  = (TOTAL | COEFF),
        };

        static const Enum<mode> modeNames;

        //- Enumeration for hydrostatic contributions
        enum hydrostaticMode : unsigned
        {
            NONE = 0,
            ADD,
            SUBTRACT
        };

        static const Enum<hydrostaticMode> hydrostaticModeNames;


private:

    // Private data

        //- Calculation mode
        mode mode_;

        //- Hydrostatic constribution mode
        hydrostaticMode hydrostaticMode_;

        //- Name of velocity field, default is "U"
        word UName_;

        //- Name of density field, default is "rho"
        word rhoName_;


        // Total pressure calculation

            //- Reference pressure level
            scalar pRef_;


        // Pressure coefficient calculation

            //- Freestream pressure
            scalar pInf_;

            //- Freestream velocity
            vector UInf_;

            //- Freestream density
            scalar rhoInf_;

            //- Flag to show whether rhoInf has been initialised
            bool rhoInfInitialised_;


        //- p +/- rgh calculation

            //- Gravity vector
            mutable dimensionedVector g_;

            //- Flag to show whether g has been initialised
            bool gInitialised_;

            //- Reference height
            mutable dimensionedScalar hRef_;

            //- Flag to show whether hRef has been initialised
            bool hRefInitialised_;


    // Private Member Functions

        //- Return the name of the derived pressure field
        word resultName() const;

        //- Multiply the static pressure p by rhoInf if necessary and return
        tmp<volScalarField> rhoScale(const volScalarField& p) const;

        //- Multiply the given field by rho or rhoInf as appropriate and return
        tmp<volScalarField> rhoScale
        (
            const volScalarField& p,
            const tmp<volScalarField>& tsf
        ) const;

        //- Add the hydrostatic contribution
        void addHydrostaticContribution
        (
            const volScalarField& p,
            volScalarField& prgh
        ) const;

        //- Calculate and return the pressure
        tmp<volScalarField> calcPressure
        (
            const volScalarField& p,
            const tmp<volScalarField>& tp
        ) const;

        //- Convert to coeff by applying the freestream dynamic pressure scaling
        tmp<volScalarField> coeff(const tmp<volScalarField>& tp) const;

        //- Calculate the derived pressure field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("pressure");


    // Constructors

        //- Construct from Time and dictionary
        pressure
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~pressure() = default;


    // Member Functions

        //- Read the pressure data
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
