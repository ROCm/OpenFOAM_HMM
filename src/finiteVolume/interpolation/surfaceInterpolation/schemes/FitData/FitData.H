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
    Foam::FitData

Description
    Data for the upwinded and centred polynomial fit interpolation schemes

SourceFiles
    FitData.C

\*---------------------------------------------------------------------------*/

#ifndef FitData_H
#define FitData_H

#include "MeshObject.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class FitData Declaration
\*---------------------------------------------------------------------------*/

template<class FitDataType, class ExtendedStencil, class Polynomial>
class FitData
:
    public MeshObject<fvMesh, FitDataType>
{
    // Private data

        //- The stencil the fit is based on
        const ExtendedStencil& stencil_;

        //- Factor the fit is allowed to deviate from linear.
        //  This limits the amount of high-order correction and increases
        //  stability on bad meshes
        const scalar linearLimitFactor_;

        //- Weights for central stencil
        const scalar centralWeight_;

        //- Dimensionality of the geometry
        const label dim_;

        //- Minimum stencil size
        const label minSize_;


    // Private member functions

        //- Find the normal direction (i) and j and k directions for face faci
        void findFaceDirs
        (
            vector& idir,        // value changed in return
            vector& jdir,        // value changed in return
            vector& kdir,        // value changed in return
            const label faci
        );

        //- Calculate the fit for the all the mesh faces
        //  and set the coefficients
        // virtual void calcFit();


public:

    //TypeName("FitData");


    // Constructors

        //- Construct from components
        FitData
        (
            const fvMesh& mesh,
            const ExtendedStencil& stencil,
            const scalar linearLimitFactor,
            const scalar centralWeight
        );


    //- Destructor
    virtual ~FitData()
    {}


    // Member functions

        //- Return reference to the stencil
        const ExtendedStencil& stencil() const
        {
            return stencil_;
        }

        //- Calculate the fit for the specified face and set the coefficients
        void calcFit
        (
            scalarList& coeffsi, // coefficients to be set
            const List<point>&,  // Stencil points
            const scalar wLin,   // Linear weight
            const label faci     // Current face index
        );

        //- Calculate the fit for all the faces
        virtual void calcFit() = 0;


        //- Delete the data when the mesh moves not implemented
        bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "FitData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
