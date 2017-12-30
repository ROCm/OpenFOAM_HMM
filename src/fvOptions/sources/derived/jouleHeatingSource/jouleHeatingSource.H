/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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
    Foam::fv::jouleHeatingSource

Group
    grpFvOptionsSources

Description
    Evolves an electrical potential equation

    \f[
        \grad \left( \sigma \grad V \right)
    \f]

    where \f$ V \f$ is electrical potential and \f$\sigma\f$ is the
    electrical current

    To provide a Joule heating contribution according to:

    Differential form of Joule heating - power per unit volume:

    \f[
        \frac{d(P)}{d(V)} = J \cdot E
    \f]

    where \f$ J \f$ is the current density and \f$ E \f$ the electric field.
    If no magnetic field is present:

    \f[
        J = \sigma E
    \f]

    The electric field given by

    \f[
        E = \grad V
    \f]

    Therefore:

    \f[
        \frac{d(P)}{d(V)} = J \cdot E
                          = (sigma E) \cdot E
                          = (sigma \grad V) \cdot \grad V
    \f]


Usage
    Isotropic (scalar) electrical conductivity
    \verbatim
    jouleHeatingSourceCoeffs
    {
        anisotropicElectricalConductivity no;

        // Optionally specify the conductivity as a function of temperature
        // Note: if not supplied, this will be read from the time directory
        sigma           table
        (
            (273        1e5)
            (1000       1e5)
        );
    }
    \endverbatim

    Anisotropic (vectorial) electrical conductivity
    jouleHeatingSourceCoeffs
    {
        anisotropicElectricalConductivity yes;

        coordinateSystem
        {
            type        cartesian;
            origin      (0 0 0);

            coordinateRotation
            {
                type        axesRotation;
                e1          (1 0 0);
                e3          (0 0 1);
            }
        }

        // Optionally specify sigma as a function of temperature
        //sigma           (31900 63800 127600);
        //
        //sigma           table
        //(
        //    (0      (0 0 0))
        //    (1000   (127600 127600 127600))
        //);
    }


    Where:
    \table
        Property     | Description               | Required  | Default value
        T            | Name of temperature field | no        | T
        sigma        | Electrical conductivity as a function of temperature |no|
        anisotropicElectricalConductivity | Anisotropic flag | yes |
    \endtable

    The electrical conductivity can be specified using either:
    - If the \c sigma entry is present the electrical conductivity is specified
      as a function of temperature using a Function1 type
    - If not present the sigma field will be read from file
    - If the anisotropicElectricalConductivity flag is set to 'true', sigma
      should be specified as a vector quantity

SourceFiles
    jouleHeatingSource.C

SeeAlso
    Foam::Function1

\*---------------------------------------------------------------------------*/

#ifndef fv_jouleHeatingSource_H
#define fv_jouleHeatingSource_H

#include "fvOption.H"
#include "Function1.H"
#include "coordinateSystem.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class jouleHeatingSource Declaration
\*---------------------------------------------------------------------------*/

class jouleHeatingSource
:
    public option
{
    // Private data

        //- Name of electrical conductivity field
        static const word sigmaName;

        //- Name of temperature field - default = "T" (optional)
        word TName_;

        //- Electrical potential field / [V]
        volScalarField V_;

        //- Flag to indicate that the electrical conductivity is anisotropic
        bool anisotropicElectricalConductivity_;

        //- Electrical conductivity as a scalar function of temperature
        autoPtr<Function1<scalar>> scalarSigmaVsTPtr_;

        //- Electrical conductivity as a vector function of temperature
        autoPtr<Function1<vector>> vectorSigmaVsTPtr_;

        //- Co-ordinate system - used for vectorial electrical conductivity
        autoPtr<coordinateSystem> coordSysPtr_;

        //- Current time index (used for updating)
        label curTimeIndex_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        jouleHeatingSource(const jouleHeatingSource&);

        //- Disallow default bitwise assignment
        void operator=(const jouleHeatingSource&);

        //- Return the co-ordinate system for anisotropic electrical
        //  conductivity
        const coordinateSystem& coordSys() const;

        //- Transform the anisotropic electrical conductivity into global system
        tmp<volSymmTensorField> transformSigma
        (
            const volVectorField& sigmaLocal
        ) const;

        //- Initialise the electrical conductivity field
        template<class Type>
        void initialiseSigma
        (
            const dictionary& dict,
            autoPtr<Function1<Type>>& sigmaVsTPtr
        );

        //- Update the electrical conductivity field
        template<class Type>
        const GeometricField<Type, fvPatchField, volMesh>&
        updateSigma(const autoPtr<Function1<Type>>& sigmaVsTPtr) const;


public:

    //- Runtime type information
    TypeName("jouleHeatingSource");


    // Constructors

        //- Construct from explicit source name and mesh
        jouleHeatingSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~jouleHeatingSource();


    // Member Functions

        // Evaluate

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "jouleHeatingSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
