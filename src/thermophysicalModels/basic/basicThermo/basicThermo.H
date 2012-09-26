/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::basicThermo

Description
    Abstract base-class for fluid and solid thermodynamic properties

SourceFiles
    basicThermo.C
    basicThermoNew.C

\*---------------------------------------------------------------------------*/

#ifndef basicThermo_H
#define basicThermo_H

#include "volFields.H"
#include "typeInfo.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "wordIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class basicThermo Declaration
\*---------------------------------------------------------------------------*/

class basicThermo
:
    public IOdictionary
{

protected:

    // Protected data

        // Fields

            //- Pressure [Pa]
            volScalarField p_;

            //- Temperature [K]
            volScalarField T_;

            //- Laminar thermal diffusuvity [kg/m/s]
            volScalarField alpha_;

        //- Should the dpdt term be included in the enthalpy equation
        Switch dpdt_;

        //- Construct as copy (not implemented)
        basicThermo(const basicThermo&);


public:

    //- Runtime type information
    TypeName("basicThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicThermo,
        fvMesh,
        (const fvMesh& mesh),
        (mesh)
    );


    // Constructors

        //- Construct from mesh
        basicThermo(const fvMesh&);

        //- Construct from mesh
        basicThermo(const fvMesh&, const dictionary&);


    //- Selector
    static autoPtr<basicThermo> New(const fvMesh&);


    //- Destructor
    virtual ~basicThermo();


    // Member functions

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate
        (
            const word& app,
            const word&
        ) const;

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate
        (
            const word& app,
            const word&,
            const word&
        ) const;

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate
        (
            const word& app,
            const word&,
            const word&,
            const word&
        ) const;

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate
        (
            const word& app,
            const word&,
            const word&,
            const word&,
            const word&
        ) const;

        //- Split name of thermo package into a list of the components names
        static wordList splitThermoName
        (
            const word& thermoName,
            const int nCmpt
        );

        //- Update properties
        virtual void correct() = 0;

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const = 0;

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const = 0;

        //- Should the dpdt term be included in the enthalpy equation
        Switch dpdt() const
        {
            return dpdt_;
        }


        // Access to thermodynamic state variables

            //- Pressure [Pa]
            //  Non-const access allowed for transport equations
            virtual volScalarField& p();

            //- Pressure [Pa]
            virtual const volScalarField& p() const;

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const = 0;

            //- Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& he() = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& he() const = 0;

            //- Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const = 0;

            //- Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Chemical enthalpy [J/kg]
            virtual tmp<volScalarField> hc() const = 0;

            //- Temperature from enthalpy/internal energy for cell-set
            virtual tmp<scalarField> THE
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const labelList& cells
            ) const = 0;

            //- Temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> THE
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const = 0;


        // Fields derived from thermodynamic state variables

            //- Temperature [K]
            virtual const volScalarField& T() const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp() const = 0;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const = 0;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- gamma = Cp/Cv []
            virtual tmp<volScalarField> gamma() const = 0;

            //- gamma = Cp/Cv for patch []
            virtual tmp<scalarField> gamma
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Heat capacity at constant pressure/volume [J/kg/K]
            virtual tmp<volScalarField> Cpv() const = 0;

            //- Heat capacity at constant pressure/volume for patch [J/kg/K]
            virtual tmp<scalarField> Cpv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Heat capacity ratio []
            virtual tmp<volScalarField> CpByCpv() const = 0;

            //- Heat capacity ratio for patch []
            virtual tmp<scalarField> CpByCpv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const = 0;


        // Access to transport state variables

            //- Thermal diffusivity for enthalpy of mixture [kg/m/s]
            virtual const volScalarField& alpha() const;


        // Fields derived from transport state variables

            //- Thermal diffusivity for temperature of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappa() const = 0;

            //- Thermal diffusivity of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappa
            (
                const label patchi
            ) const = 0;

            //- Effective thermal diffusivity of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappaEff
            (
                const volScalarField&
            ) const = 0;

            //- Effective thermal diffusivity of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const = 0;

            //- Effective thermal diffusivity of mixture [J/m/s/K]
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const = 0;

            //- Effective thermal diffusivity of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const = 0;


        //- Read thermophysicalProperties dictionary
        virtual bool read();
};


template<class Thermo>
autoPtr<Thermo> NewThermo
(
    const fvMesh& mesh
)
{
    IOdictionary thermoDict
    (
        IOobject
        (
            "thermophysicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    word thermoTypeName;

    if (thermoDict.isDict("thermoType"))
    {
        const dictionary& thermoTypeDict(thermoDict.subDict("thermoType"));

        Info<< "Selecting thermodynamics package " << thermoTypeDict << endl;

        const int nCmpt = 7;
        const char* cmptNames[nCmpt] =
        {
            "type",
            "mixture",
            "transport",
            "thermo",
            "equationOfState",
            "specie",
            "energy"
        };

        // Construct the name of the thermo package from the components
        thermoTypeName =
            word(thermoTypeDict.lookup("type")) + '<'
          + word(thermoTypeDict.lookup("mixture")) + '<'
          + word(thermoTypeDict.lookup("transport")) + '<'
          + word(thermoTypeDict.lookup("thermo")) + '<'
          + word(thermoTypeDict.lookup("equationOfState")) + '<'
          + word(thermoTypeDict.lookup("specie")) + ">>,"
          + word(thermoTypeDict.lookup("energy")) + ">>>";

        // Lookup the thermo package
        typename Thermo::fvMeshConstructorTable::iterator cstrIter =
            Thermo::fvMeshConstructorTablePtr_->find(thermoTypeName);

        // Print error message if package not found in the table
        if (cstrIter == Thermo::fvMeshConstructorTablePtr_->end())
        {
            FatalErrorIn(Thermo::typeName + "::New(const fvMesh&)")
                << "Unknown " << Thermo::typeName << " type " << nl
                << "thermoType" << thermoTypeDict << nl << nl
                << "Valid " << Thermo::typeName << " types are:" << nl << nl;

            // Get the list of all the suitable thermo packages available
            wordList validThermoTypeNames
            (
                Thermo::fvMeshConstructorTablePtr_->sortedToc()
            );

            // Build a table of the thermo packages constituent parts
            // Note: row-0 contains the names of constituent parts
            List<wordList> validThermoTypeNameCmpts
            (
                validThermoTypeNames.size() + 1
            );

            validThermoTypeNameCmpts[0].setSize(nCmpt);
            forAll(validThermoTypeNameCmpts[0], j)
            {
                validThermoTypeNameCmpts[0][j] = cmptNames[j];
            }

            // Split the thermo package names into their constituent parts
            forAll(validThermoTypeNames, i)
            {
                validThermoTypeNameCmpts[i+1] =
                    Thermo::splitThermoName(validThermoTypeNames[i], nCmpt);
            }

            // Print the table of available packages
            // in terms of their constituent parts
            printTable(validThermoTypeNameCmpts, FatalError);

            FatalError<< exit(FatalError);
        }

        return autoPtr<Thermo>(cstrIter()(mesh));
    }
    else
    {
        thermoTypeName = word(thermoDict.lookup("thermoType"));

        Info<< "Selecting thermodynamics package " << thermoTypeName << endl;

        typename Thermo::fvMeshConstructorTable::iterator cstrIter =
            Thermo::fvMeshConstructorTablePtr_->find(thermoTypeName);

        if (cstrIter == Thermo::fvMeshConstructorTablePtr_->end())
        {
            FatalErrorIn(Thermo::typeName + "::New(const fvMesh&)")
                << "Unknown " << Thermo::typeName << " type "
                << thermoTypeName << nl << nl
                << "Valid " << Thermo::typeName << " types are:" << nl
                << Thermo::fvMeshConstructorTablePtr_->sortedToc() << nl
                << exit(FatalError);
        }

        return autoPtr<Thermo>(cstrIter()(mesh));
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
