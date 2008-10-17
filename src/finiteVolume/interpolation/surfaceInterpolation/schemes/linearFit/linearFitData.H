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
    linearFitData

Description
    Data for the linear fit correction interpolation scheme

SourceFiles
    linearFitData.C

\*---------------------------------------------------------------------------*/

#ifndef linearFitData_H
#define linearFitData_H

#include "MeshObject.H"
#include "fvMesh.H"
#include "extendedStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class globalIndex;

/*---------------------------------------------------------------------------*\
                    Class linearFitData Declaration
\*---------------------------------------------------------------------------*/

class linearFitData
:
    public MeshObject<fvMesh, linearFitData>
{
    // Private data

        //- weights for central stencil
        const scalar centralWeight_;

        //- dimensionality of the geometry
        const label dim_;

        //- minimum stencil size
        const label minSize_;

        //- Extended stencil addressing
        extendedStencil stencil_;

        //- For each cell in the mesh store the values which multiply the
        //  values of the stencil to obtain the gradient for each direction
        List<scalarList> fit_;


    // Private member functions

        //- Find the normal direction and i, j and k directions for face faci
        static void findFaceDirs
        (
            vector& idir,        // value changed in return
            vector& jdir,        // value changed in return
            vector& kdir,        // value changed in return
            const fvMesh& mesh,
            const label faci
        );

        label calcFit(const List<point>&, const label faci);


public:

    TypeName("linearFitData");


    // Constructors

        explicit linearFitData
        (
            const fvMesh& mesh,
            scalar cWeightDim
        );


    // Destructor

        virtual ~linearFitData()
        {}


    // Member functions


        //- Return reference to the stencil
        const extendedStencil& stencil() const
        {
            return stencil_;
        }

        //- Return reference to fit coefficients
        const List<scalarList>& fit() const
        {
            return fit_;
        }

        //- Delete the data when the mesh moves not implemented
        virtual bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
