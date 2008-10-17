/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::quadraticFit

Description
    Quadratic fit interpolation scheme which applies an explicit correction to
    linear.

SourceFiles
    quadraticFit.C

\*---------------------------------------------------------------------------*/

#ifndef quadraticFit_H
#define quadraticFit_H

#include "quadraticFitData.H"
#include "linear.H"
#include "centredCFCStencilObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class quadraticFit Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class quadraticFit
:
    public linear<Type>
{
    // Private Data

        //- Factor the fit is allowed to deviate from linear.
        //  This limits the amount of high-order correction and increases
        //  stability on bad meshes
        const scalar linearLimitFactor_;

        //- Weights for central stencil
        const scalar centralWeight_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        quadraticFit(const quadraticFit&);

        //- Disallow default bitwise assignment
        void operator=(const quadraticFit&);


public:

    //- Runtime type information
    TypeName("quadraticFit");


    // Constructors

        //- Construct from mesh and Istream
        quadraticFit(const fvMesh& mesh, Istream& is)
        :
            linear<Type>(mesh),
            linearLimitFactor_(readScalar(is)),
            centralWeight_(readScalar(is))
        {}


        //- Construct from mesh, faceFlux and Istream
        quadraticFit
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            linear<Type>(mesh),
            linearLimitFactor_(readScalar(is)),
            centralWeight_(readScalar(is))
        {}


    // Member Functions

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
            const fvMesh& mesh = this->mesh();

            const extendedCentredStencil& stencil =
            centredCFCStencilObject::New
            (
                mesh
            );

            const quadraticFitData& cfd = quadraticFitData::New
            (
                mesh,
                stencil,
                linearLimitFactor_,
                centralWeight_
            );

            const List<scalarList>& f = cfd.fit();

            return stencil.weightedSum(vf, f);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
