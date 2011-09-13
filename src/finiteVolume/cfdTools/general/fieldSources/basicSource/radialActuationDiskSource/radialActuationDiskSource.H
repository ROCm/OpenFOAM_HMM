/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::radialActuationDiskSource

Description
    Actuation disk zone definition.
    Constant values for momentum source for actuation disk

    T = 2*rho*A*sqr(Uo)*a*(1-a)
    U1 = (1 -a)Uo
    where:
          A: disk area
          Uo: upstream velocity
          a: 1 - Cp/Ct
          U1: velocity at the disk


    The thrust is distributed by a radial function:
    thrust(r) = T*(C0 + C1*r^2 + C2*r^4)


SourceFiles
    radialActuationDiskSource.C
    radialActuationDiskSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef radialActuationDiskSource_H
#define radialActuationDiskSource_H

#include "DimensionedField.H"
#include "volFieldsFwd.H"
#include "FixedList.H"
#include "actuationDiskSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class radialActuationDiskSource Declaration
\*---------------------------------------------------------------------------*/

class radialActuationDiskSource
:
    public actuationDiskSource
{
    // Private data

        //- Sub dictionary with model information
        const dictionary& dict_;

        //- Coeffcients for the radial distribution
        FixedList<scalar, 3> coeffs_;


    // Private Member Functions

        //- Add resistance to the UEqn
        template<class RhoFieldType>
        void addRadialActuationDiskAxialInertialResistance
        (
            vectorField& Usource,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;

        //- Disallow default bitwise copy construct
        radialActuationDiskSource(const radialActuationDiskSource&);

        //- Disallow default bitwise assignment
        void operator=(const radialActuationDiskSource&);


public:

    //- Runtime type information
    TypeName("radialActuationDiskSource");


    // Constructors

        //- Construct from components
        radialActuationDiskSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~radialActuationDiskSource()
    {}


        // Public Functions

            //-Source term to fvMatrix<vector>
            virtual void addSu(fvMatrix<vector>& UEqn);


        // I-O

            //- Write data
            virtual void writeData(Ostream&) const;

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "radialActuationDiskSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
