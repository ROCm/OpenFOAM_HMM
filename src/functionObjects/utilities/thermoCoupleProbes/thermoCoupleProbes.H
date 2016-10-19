/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD.Ltd
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
    Foam::functionObjects::thermoCoupleProbes

Description
    Sample probe for temperature using a thermocouple. Uses the correlation:

    Nu = 2.0 + (0.4*sqrt(Re) + 0.06*pow(Re, 2/3))*pow(Pr, 0.4)

    Example of function object specification:
    \verbatim
    probes
    {
        type            thermoCoupleProbes;

        solver          rodas23;
        absTol          1e-12;
        relTol          1e-8;

        interpolationScheme cellPoint;

        // thermocouple properties
        rho             8908;
        Cp              440;
        d               1e-3;
        epsilon         0.85;

        radName         G;

        functionObjectLibs ("libsampling.so");
        outputControl   timeStep;
        outputInterval  1;
        probeLocations
        (
            ( 0.5 0.5 0.5 )
        );
        fields
        (
            T
        );
    }
    \endverbatim


SourceFiles
    thermoCoupleProbes.C

\*---------------------------------------------------------------------------*/

#ifndef thermoCoupleProbes_H
#define thermoCoupleProbes_H

#include "probes.H"
#include "ODESystem.H"
#include "ODESolver.H"
#include "basicThermo.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class thermoCoupleProbes Declaration
\*---------------------------------------------------------------------------*/

class thermoCoupleProbes
:
    public probes,
    public ODESystem
{
protected:

    // Protected data

        //- Thermocouple density
        scalar rho_;

        //- Thermocouple heat capacity
        scalar Cp_;

        //- Thermocouple diameter
        scalar d_;

        //- Thermocouple emissivity
        scalar epsilon_;

        //- Name for the incident radiation field
        word radName_;

         //- Fluid thermo reference
        const fluidThermo& thermo_;

        //- ODESolver
        autoPtr<ODESolver> odeSolver_;

        //- Cached thermocouple temperature
        scalarField Ttc_;


    // Protected Member Functions

        //- Sample and write a particular volume field
        template<class Type>
        void sampleAndWrite
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );


        //- Read dictionary settings
        void readDict(const dictionary& dict);


private:

        //- Disallow default bitwise copy construct
        thermoCoupleProbes(const thermoCoupleProbes&);

        //- Disallow default bitwise assignment
        void operator=(const thermoCoupleProbes&);


public:

    //- Runtime type information
    TypeName("thermoCoupleProbes");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        thermoCoupleProbes
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const bool loadFromFiles = false,
            const bool findElements = true
        );


    //- Destructor
    virtual ~thermoCoupleProbes();


    // ODE functions (overriding abstract functions in ODE.H)

    //- Number of ODE's to solve
    virtual label nEqns() const;

    virtual void derivatives
    (
        const scalar t,
        const scalarField& c,
        scalarField& dcdt
    ) const;

    virtual void jacobian
    (
        const scalar t,
        const scalarField& c,
        scalarField& dcdt,
        scalarSquareMatrix& dfdc
    ) const;

    //- Public members

        //- Sample and write
        virtual bool write();

        //- Execute, currently does nothing
        virtual bool execute();

        //- Read
        virtual bool read(const dictionary&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "thermoCoupleProbesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
