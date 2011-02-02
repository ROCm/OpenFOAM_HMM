/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
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
    quadraticUpwind

Description
    quadraticUpwind interpolation scheme class derived from linearUpwind and
    returns blended linear/upwind weighting factors and also applies a explicit
    gradient-based correction obtained from the linearUpwind scheme.

SourceFiles
    quadraticUpwind.C

\*---------------------------------------------------------------------------*/

#ifndef quadraticUpwind_H
#define quadraticUpwind_H

#include "linearUpwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class quadraticUpwind Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class quadraticUpwind
:
    public linearUpwind<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        quadraticUpwind(const quadraticUpwind&);

        //- Disallow default bitwise assignment
        void operator=(const quadraticUpwind&);


public:

    //- Runtime type information
    TypeName("quadraticUpwind");


    // Constructors

        //- Construct from mesh and Istream
        quadraticUpwind
        (
            const fvMesh& mesh,
            Istream& schemeData
        )
        :
            linearUpwind<Type>(mesh, schemeData)
        {}

        //- Construct from mesh, faceFlux and Istream
        quadraticUpwind
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        )
        :
            linearUpwind<Type>(mesh, faceFlux, schemeData)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            return
                0.75*this->mesh().surfaceInterpolation::weights()
              + 0.25*linearUpwind<Type>::weights();
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            return 0.25*linearUpwind<Type>::correction(vf);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
