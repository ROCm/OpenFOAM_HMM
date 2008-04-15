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

#include "channelIndex.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

channelIndex::channelIndex(const fvMesh& m)
:
    indexingDict_
    (
        IOobject
        (
            "postChannelDict",
            m.time().constant(),
            m,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    nx_(readLabel(indexingDict_.lookup("Nx"))),
    ny_(indexingDict_.lookup("Ny")),
    nz_(readLabel(indexingDict_.lookup("Nz"))),
    symmetric_
    (
        readBool(indexingDict_.lookup("symmetric"))
    ),
    cumNy_(ny_.size()),
    nLayers_(ny_[0])
{
    // initialise the layers
    cumNy_[0] = ny_[0];

    for (label j=1; j<ny_.size(); j++)
    {
        nLayers_ += ny_[j];
        cumNy_[j] = ny_[j]+cumNy_[j-1];
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

channelIndex::~channelIndex()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalarField channelIndex::collapse
(
    const volScalarField& vsf,
    const bool asymmetric
) const
{
    scalarField cs(nLayers(), 0.0);

    forAll(cs, j)
    {
        // sweep over all cells in this layer
        for (label i=0; i<nx(); i++)
        {
            for (label k=0; k<nz(); k++)
            {
                cs[j] += vsf[operator()(i,j,k)];
            }
        }

        // and divide by the number of cells in the layer
        cs[j] /= scalar(nx()*nz());
    }

    if (symmetric_)
    {
        label nlb2 = nLayers()/2;

        if (asymmetric)
        {
            for (label j=0; j<nlb2; j++)
            {
                cs[j] = 0.5*(cs[j] - cs[nLayers() - j - 1]);
            }
        }
        else
        {
            for (label j=0; j<nlb2; j++)
            {
                cs[j] = 0.5*(cs[j] + cs[nLayers() - j - 1]);
            }
        }

        cs.setSize(nlb2);
    }

    return cs;
}


scalarField channelIndex::y
(
    const volVectorField& cellCentres
) const
{
    if (symmetric_)
    {
        scalarField Y(nLayers()/2);

        for (label j=0; j<nLayers()/2; j++)
        {
            Y[j] = cellCentres[operator()(0, j, 0)].y();
        }

        return Y;
    }
    else
    {
        scalarField Y(nLayers());

        for (label j=0; j<nLayers(); j++)
        {
            Y[j] = cellCentres[operator()(0, j, 0)].y();
        }

        return Y;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

label channelIndex::operator()
(
    const label Jx,
    const label Jy,
    const label Jz
) const
{
    label index(0);

    // count up `full' layers in the mesh
    label j(0);
    label tmpJy(Jy);

    while(Jy >= cumNy_[j])
    {
        index += nx_*ny_[j]*nz_;
        tmpJy -= ny_[j];
        j++;
    }

    index += Jx + nx_*tmpJy + nx_*ny_[j]*Jz;

    return index;
}


// ************************************************************************* //
