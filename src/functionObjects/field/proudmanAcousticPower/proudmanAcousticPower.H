/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::functionObjects::proudmanAcousticPower

Group
    grpFieldFunctionObjects

Description
    Calculates the acoustic power due to the volume of isotropic turbulence
    using Proudman's formula

    The acoustic power \f$ P_A \f$ [W/m3] in terms of turbulence \f$ k \f$
    and \f$ \epsilon \f$ is given as:

        \f[
            P_A = alpha_\epsilon \rho \epsilon M_t^5
        \f]

    where \f$ alpha_\epsilon \f$ is a constant (0.1) and

        \f[
            M_t = \frac{\sqrt{2 k}}{a_0}
        \f]

    with \f$ a_0 \f$ the speed of sound.  The acoustic power is also output in
    dB using:

        \f[
            L_P = 10 \log \frac{P_A}{P_ref}
        \f]

    where \f$ P_ref \f$ is a constant (1e-12 W/m3)

Usage
    Example of function object specification to calculate the Proudman acoustic
    power:
    \verbatim
    proudmanAcousticPower1
    {
        type        proudmanAcousticPower;
        libs        ("libfieldFunctionObjects.so");
        ...

        // Required additional entries for incompressible calculations
        rhoInf      1.225;
        aRef        340;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property | Description                              | Required | Default
        type     | type name: proudmanAcousticPower         | yes |
        rhoInf   | Freestream density (for incompressible)  | no  |
        aRef     | Reference speed of sound (for incompressible) | no  |
        alphaEps | Model coefficient                        | no  | 0.1
    \endtable

Note
    The freestream density and reference speed of sound are only necessary
    when a thermodynamics package is unavailable, typically for incompressible
    cases.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    proudmanAcousticPower.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_proudmanAcousticPower_H
#define functionObjects_proudmanAcousticPower_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class proudmanAcousticPower Declaration
\*---------------------------------------------------------------------------*/

class proudmanAcousticPower
:
    public fvMeshFunctionObject
{
private:

    // Private Data

        //- Freestream density (incompressible calcs only)
        dimensionedScalar rhoInf_;

        //- Reference speed of sound (incompressible calcs only)
        dimensionedScalar aRef_;

        //- Model coefficient; default = 0.1
        scalar alphaEps_;


    // Private Member Functions

        //- Multiply the field by density and return
        tmp<volScalarField> rhoScale(const tmp<volScalarField>& fld) const;

        //- Speed of sound
        tmp<volScalarField> a() const;


public:

    //- Runtime type information
    TypeName("proudmanAcousticPower");


    // Constructors

        //- Construct from Time and dictionary
        proudmanAcousticPower
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~proudmanAcousticPower() = default;


    // Member Functions

        //- Read the Proudman acoustic power data
        virtual bool read(const dictionary&);

        //- Calculate the Proudman acoustic power
        virtual bool execute();

        //- Write the Proudman acoustic power
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
