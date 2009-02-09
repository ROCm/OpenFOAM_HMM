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

\*---------------------------------------------------------------------------*/

#include "UpwindFitData.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "SVD.H"
#include "syncTools.H"
#include "extendedUpwindStencil.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Polynomial>
Foam::UpwindFitData<Polynomial>::UpwindFitData
(
    const fvMesh& mesh,
    const extendedUpwindStencil& stencil,
    const scalar linearLimitFactor,
    const scalar centralWeight
)
:
    FitData<UpwindFitData<Polynomial>, extendedUpwindStencil, Polynomial>
    (
        mesh, stencil, linearLimitFactor, centralWeight
    ),
    owncoeffs_(mesh.nFaces()),
    neicoeffs_(mesh.nFaces())
{
    if (debug)
    {
        Info<< "Contructing UpwindFitData<Polynomial>" << endl;
    }

    calcFit();

    if (debug)
    {
        Info<< "UpwindFitData<Polynomial>::UpwindFitData() :"
            << "Finished constructing polynomialFit data"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Polynomial>
void Foam::UpwindFitData<Polynomial>::calcFit()
{
    const fvMesh& mesh = this->mesh();

    // Get the cell/face centres in stencil order.
    // Upwind face stencils no good for triangles or tets.
    // Need bigger stencils
    List<List<point> > ownStencilPoints(mesh.nFaces());
    this->stencil().collectData
    (
        this->stencil().ownMap(),
        this->stencil().ownStencil(),
        mesh.C(),
        ownStencilPoints
    );
    List<List<point> > neiStencilPoints(mesh.nFaces());
    this->stencil().collectData
    (
        this->stencil().neiMap(),
        this->stencil().neiStencil(),
        mesh.C(),
        neiStencilPoints
    );

    // find the fit coefficients for every owner and neighbour of ever face

    const surfaceScalarField& w = mesh.surfaceInterpolation::weights();

    for(label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        FitData<UpwindFitData<Polynomial>, extendedUpwindStencil, Polynomial>::
           calcFit(owncoeffs_[facei], ownStencilPoints[facei], w[facei], facei);
        FitData<UpwindFitData<Polynomial>, extendedUpwindStencil, Polynomial>::
           calcFit(neicoeffs_[facei], neiStencilPoints[facei], w[facei], facei);
    }

    const surfaceScalarField::GeometricBoundaryField& bw = w.boundaryField();

    forAll(bw, patchi)
    {
        const fvsPatchScalarField& pw = bw[patchi];

        if (pw.coupled())
        {
            label facei = pw.patch().patch().start();

            forAll(pw, i)
            {
                FitData
                <
                    UpwindFitData<Polynomial>, extendedUpwindStencil, Polynomial
                >::calcFit
                (
                    owncoeffs_[facei], ownStencilPoints[facei], pw[i], facei
                );
                FitData
                <
                    UpwindFitData<Polynomial>, extendedUpwindStencil, Polynomial
                >::calcFit
                (
                    neicoeffs_[facei], neiStencilPoints[facei], pw[i], facei
                );
                facei++;
            }
        }
    }
}


// ************************************************************************* //
