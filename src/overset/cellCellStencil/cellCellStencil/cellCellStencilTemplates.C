/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::cellCellStencil::interpolate
(
    Field<T>& psi,
    const fvMesh& mesh,
    const cellCellStencil& overlap,
    const List<scalarList>& wghts
)
{
    const labelListList& stencil = overlap.cellStencil();

    if (stencil.size() != mesh.nCells())
    {
        return;
    }

    const mapDistribute& map = overlap.cellInterpolationMap();
    const labelList& cellIDs = overlap.interpolationCells();
    const scalarList& factor = overlap.cellInterpolationWeight();

    Field<T> work(psi);
    map.mapDistributeBase::distribute(work, UPstream::msgType()+1);

    forAll(cellIDs, i)
    {
        const label celli = cellIDs[i];

        const scalarList& w = wghts[celli];
        const labelList& nbrs = stencil[celli];
        const scalar f = factor[celli];

        if (nbrs.size() == 0 && f != 0.0)
        {
            FatalErrorInFunction << "problem: cell:" << celli
                << " at:" << mesh.cellCentres()[celli]
                << " type:" << overlap.cellTypes()[celli]
                << " stencil:" << nbrs
                << " factor:" << f << exit(FatalError);
        }

        T s(pTraits<T>::zero);
        forAll(nbrs, nbrI)
        {
            s += w[nbrI]*work[nbrs[nbrI]];
        }

        psi[celli] = (1.0-f)*psi[celli] + f*s;
    }
}


template<class T>
void Foam::cellCellStencil::interpolate(const fvMesh& mesh, Field<T>& psi) const
{
    const cellCellStencil& overlap = *this;

    interpolate
    (
        psi,
        mesh,
        overlap,
        overlap.cellInterpolationWeights()
    );
}


template<class GeoField>
void Foam::cellCellStencil::interpolate(GeoField& psi) const
{
    interpolate
    (
        psi.mesh(),
        psi.primitiveFieldRef()
    );
    psi.correctBoundaryConditions();
}


template<class GeoField>
void Foam::cellCellStencil::interpolate
(
    const fvMesh& mesh,
    const wordHashSet& suppressed
) const
{
    const cellCellStencil& overlap = *this;

    auto flds(mesh.lookupClass<GeoField>());
    for (auto fldPtr : flds)
    {
        const word& name = fldPtr->name();
        if (!suppressed.found(baseName(name)))
        {
            if (debug)
            {
                Pout<< "cellCellStencil::interpolate: interpolating : "
                    << name << endl;
            }

            auto& fld = const_cast<GeoField&>(*fldPtr);

            interpolate
            (
                fld.primitiveFieldRef(),
                mesh,
                overlap,
                overlap.cellInterpolationWeights()
            );
        }
        else
        {
            if (debug)
            {
                Pout<< "cellCellStencil::interpolate: skipping : " << name
                    << endl;
            }
        }
    }
}

template<class Type>
Foam::tmp<Foam::volScalarField>
Foam::cellCellStencil::createField
(
    const fvMesh& mesh,
    const word& name,
    const UList<Type>& psi
)
{
    auto tfld = tmp<volScalarField>::New
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    );
    volScalarField& fld = tfld.ref();

    forAll(psi, cellI)
    {
        fld[cellI] = psi[cellI];
    }
    return tfld;
}


template<class GeoField, class SuppressBC>
void Foam::cellCellStencil::correctBoundaryConditions
(
    GeoField& psi
)
{
    // Version of GeoField::correctBoundaryConditions that exclude evaluation
    // of oversetFvPatchFields

    typename GeoField::Boundary& bfld = psi.boundaryFieldRef();

    const label nReq = Pstream::nRequests();

    forAll(bfld, patchi)
    {
        if (!isA<SuppressBC>(bfld[patchi]))
        {
            bfld[patchi].initEvaluate(Pstream::commsTypes::nonBlocking);
        }
    }

    // Block for any outstanding requests
    if (Pstream::parRun())
    {
        Pstream::waitRequests(nReq);
    }

    forAll(bfld, patchi)
    {
        if (!isA<SuppressBC>(bfld[patchi]))
        {
            bfld[patchi].evaluate(Pstream::commsTypes::nonBlocking);
        }
    }
}


// ************************************************************************* //
