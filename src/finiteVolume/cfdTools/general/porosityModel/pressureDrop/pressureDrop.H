/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::pressureDrop

Description
    Pressure drop porosity model, whereby the user speciefies the pressure
    drop per unit length as a function of mass flow rate, over a fixed distance.

    Example usage:
    \verbatim
    pressureDropCoeffs
    {
        coordinateSystem
        {
            e1  (1 0 0);
            e2  (0 1 0);
        }

        lRef    0.5;

        mDotvsDp table
        (
            (0   0)
            (0.2 20)
            (0.4 40)
            (0.6 65)
            (0.8 95)
            (1.0 125)
            (1.5 220)
            (2.0 320)
            (2.5 435)
            (3.0 560)
            (3.5 700)
        );
    }
    \endverbatim

Note
    The mDotvsDp entry is a DataEntry type, able to describe the pressure drop
    per unit length as a function mass flow rate.  The example above gives the
    usage for supplying in-line tabulated data.

SeeAlso
    Foam::DataEntry

SourceFiles
    pressureDrop.C

\*---------------------------------------------------------------------------*/

#ifndef pressureDrop_H
#define pressureDrop_H

#include "porosityModel.H"
#include "coordinateSystem.H"
#include "autoPtr.H"
#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace porosityModels
{

/*---------------------------------------------------------------------------*\
                      Class pressureDrop Declaration
\*---------------------------------------------------------------------------*/

class pressureDrop
:
    public porosityModel
{
private:

    // Private data

        //- Local co-ordinate system
        coordinateSystem coordSys_;

        //- Pressure drop per unit length [Pa/m] as a function of
        //  mass flow rate [kg/s]
        autoPtr<DataEntry<scalar> > mDotvsDp_;

        //- Distance over which pressure drop is applied
        scalar lRef_;

        //- Name of density field
        word rhoName_;


    // Private Member Functions

        //- Apply
        template<class RhoFieldType>
        void apply
        (
            scalarField& Udiag,
            vectorField& Usource,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U,
            const scalar rhoScale
        ) const;

        //- Apply
        template<class RhoFieldType>
        void apply
        (
            tensorField& AU,
            const RhoFieldType& rho,
            const vectorField& U,
            const scalar rhoScale
        ) const;

        //- Disallow default bitwise copy construct
        pressureDrop(const pressureDrop&);

        //- Disallow default bitwise assignment
        void operator=(const pressureDrop&);


public:

    //- Runtime type information
    TypeName("pressureDrop");

    //- Constructor
    pressureDrop
    (
        const word& name,
        const word& modelType,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& cellZoneName
    );

    //- Destructor
    virtual ~pressureDrop();


    // Member Functions

        //- Add resistance
        virtual void correct(fvVectorMatrix& UEqn) const;

        //- Add resistance
        virtual void correct
        (
            fvVectorMatrix& UEqn,
            const volScalarField& rho,
            const volScalarField& mu
        ) const;

        //- Add resistance
        virtual void correct
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU
        ) const;


    // I-O

        //- Write
        void writeData(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace porosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pressureDropTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
