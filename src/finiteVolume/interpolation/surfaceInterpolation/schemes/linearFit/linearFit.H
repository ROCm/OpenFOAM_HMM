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
    linearFit

Description
    Linear fit interpolation scheme which applies an explicit correction to
    linear.

SourceFiles
    linearFit.C

\*---------------------------------------------------------------------------*/

#ifndef linearFit_H
#define linearFit_H

#include "linear.H"
#include "linearFitData.H"
#include "extendedStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class linearFit Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class linearFit
:
    public linear<Type>
{
    // Private Data
        const scalar centralWeight_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        linearFit(const linearFit&);

        //- Disallow default bitwise assignment
        void operator=(const linearFit&);


public:

    //- Runtime type information
    TypeName("linearFit");


    // Constructors

        //- Construct from mesh and Istream
        linearFit(const fvMesh& mesh, Istream& is)
        :
            linear<Type>(mesh),
            centralWeight_(readScalar(is))
        {}


        //- Construct from mesh, faceFlux and Istream
        linearFit
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            linear<Type>(mesh),
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

            const linearFitData& cfd = linearFitData::New
            (
                mesh,
                centralWeight_
            );

            const extendedStencil& stencil = cfd.stencil();
            const List<scalarList>& f = cfd.fit();

            return stencil.weightedSum(vf, f);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
