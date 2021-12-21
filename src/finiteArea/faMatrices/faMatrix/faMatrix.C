/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "areaFields.H"
#include "edgeFields.H"
#include "calculatedFaPatchFields.H"
#include "zeroGradientFaPatchFields.H"
#include "UIndirectList.H"
#include "UniformList.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
template<class Type2>
void Foam::faMatrix<Type>::addToInternalField
(
    const labelUList& addr,
    const Field<Type2>& pf,
    Field<Type2>& intf
) const
{
    if (addr.size() != pf.size())
    {
        FatalErrorInFunction
            << "addressing (" << addr.size()
            << ") and field (" << pf.size() << ") are different sizes" << endl
            << abort(FatalError);
    }

    forAll(addr, faceI)
    {
        intf[addr[faceI]] += pf[faceI];
    }
}


template<class Type>
template<class Type2>
void Foam::faMatrix<Type>::addToInternalField
(
    const labelUList& addr,
    const tmp<Field<Type2>>& tpf,
    Field<Type2>& intf
) const
{
    addToInternalField(addr, tpf(), intf);
    tpf.clear();
}


template<class Type>
template<class Type2>
void Foam::faMatrix<Type>::subtractFromInternalField
(
    const labelUList& addr,
    const Field<Type2>& pf,
    Field<Type2>& intf
) const
{
    if (addr.size() != pf.size())
    {
        FatalErrorInFunction
            << "addressing (" << addr.size()
            << ") and field (" << pf.size() << ") are different sizes" << endl
            << abort(FatalError);
    }

    forAll(addr, faceI)
    {
        intf[addr[faceI]] -= pf[faceI];
    }
}


template<class Type>
template<class Type2>
void Foam::faMatrix<Type>::subtractFromInternalField
(
    const labelUList& addr,
    const tmp<Field<Type2>>& tpf,
    Field<Type2>& intf
) const
{
    subtractFromInternalField(addr, tpf(), intf);
    tpf.clear();
}


template<class Type>
void Foam::faMatrix<Type>::addBoundaryDiag
(
    scalarField& diag,
    const direction solveCmpt
) const
{
    forAll(internalCoeffs_, patchI)
    {
        addToInternalField
        (
            lduAddr().patchAddr(patchI),
            internalCoeffs_[patchI].component(solveCmpt),
            diag
        );
    }
}


template<class Type>
void Foam::faMatrix<Type>::addCmptAvBoundaryDiag(scalarField& diag) const
{
    forAll(internalCoeffs_, patchI)
    {
        addToInternalField
        (
            lduAddr().patchAddr(patchI),
            cmptAv(internalCoeffs_[patchI]),
            diag
        );
    }
}


template<class Type>
void Foam::faMatrix<Type>::addBoundarySource
(
    Field<Type>& source,
    const bool couples
) const
{
    forAll(psi_.boundaryField(), patchI)
    {
        const faPatchField<Type>& ptf = psi_.boundaryField()[patchI];
        const Field<Type>& pbc = boundaryCoeffs_[patchI];

        if (!ptf.coupled())
        {
            addToInternalField(lduAddr().patchAddr(patchI), pbc, source);
        }
        else if (couples)
        {
            tmp<Field<Type>> tpnf = ptf.patchNeighbourField();
            const Field<Type>& pnf = tpnf();

            const labelUList& addr = lduAddr().patchAddr(patchI);

            forAll(addr, facei)
            {
                source[addr[facei]] += cmptMultiply(pbc[facei], pnf[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template<class Type>
Foam::faMatrix<Type>::faMatrix
(
    const GeometricField<Type, faPatchField, areaMesh>& psi,
    const dimensionSet& ds
)
:
    lduMatrix(psi.mesh()),
    psi_(psi),
    dimensions_(ds),
    source_(psi.size(), Zero),
    internalCoeffs_(psi.mesh().boundary().size()),
    boundaryCoeffs_(psi.mesh().boundary().size()),
    faceFluxCorrectionPtr_(nullptr)
{
    DebugInFunction
        << "constructing faMatrix<Type> for field " << psi_.name()
        << endl;

    // Initialise coupling coefficients
    forAll(psi.mesh().boundary(), patchI)
    {
        internalCoeffs_.set
        (
            patchI,
            new Field<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                Zero
            )
        );

        boundaryCoeffs_.set
        (
            patchI,
            new Field<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                Zero
            )
        );
    }

    // Update the boundary coefficients of psi without changing its event No.
    auto& psiRef =
        const_cast<GeometricField<Type, faPatchField, areaMesh>&>(psi_);

    const label currentStatePsi = psiRef.eventNo();
    psiRef.boundaryFieldRef().updateCoeffs();
    psiRef.eventNo() = currentStatePsi;
}


template<class Type>
Foam::faMatrix<Type>::faMatrix(const faMatrix<Type>& fam)
:
    refCount(),
    lduMatrix(fam),
    psi_(fam.psi_),
    dimensions_(fam.dimensions_),
    source_(fam.source_),
    internalCoeffs_(fam.internalCoeffs_),
    boundaryCoeffs_(fam.boundaryCoeffs_),
    faceFluxCorrectionPtr_(nullptr)
{
    DebugInFunction
        << "copying faMatrix<Type> for field " << psi_.name()
        << endl;

    if (fam.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = new
        GeometricField<Type, faePatchField, edgeMesh>
        (
            *(fam.faceFluxCorrectionPtr_)
        );
    }
}


template<class Type>
Foam::faMatrix<Type>::faMatrix
(
    const GeometricField<Type, faPatchField, areaMesh>& psi,
    Istream& is
)
:
    lduMatrix(psi.mesh()),
    psi_(psi),
    dimensions_(is),
    source_(is),
    internalCoeffs_(psi.mesh().boundary().size()),
    boundaryCoeffs_(psi.mesh().boundary().size()),
    faceFluxCorrectionPtr_(nullptr)
{
    DebugInFunction
        << "constructing faMatrix<Type> for field " << psi_.name()
        << endl;

    // Initialise coupling coefficients
    forAll(psi.mesh().boundary(), patchI)
    {
        internalCoeffs_.set
        (
            patchI,
            new Field<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                Zero
            )
        );

        boundaryCoeffs_.set
        (
            patchI,
            new Field<Type>
            (
                psi.mesh().boundary()[patchI].size(),
                Zero
            )
        );
    }
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::faMatrix<Type>::clone() const
{
    return tmp<faMatrix<Type>>::New(*this);
}


// * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //


template<class Type>
Foam::faMatrix<Type>::~faMatrix()
{
    DebugInFunction
        << "Destroying faMatrix<Type> for field " << psi_.name() << endl;

    deleteDemandDrivenData(faceFluxCorrectionPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
template<template<class> class ListType>
void Foam::faMatrix<Type>::setValuesFromList
(
    const labelUList& faceLabels,
    const ListType<Type>& values
)
{
#if 0  /* Specific to foam-extend */
    // Record face labels of eliminated equations
    for (const label i : faceLabels)
    {
        this->eliminatedEqns().insert(i);
    }
#endif

    const faMesh& mesh = psi_.mesh();

    const labelListList& edges = mesh.patch().faceEdges();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    scalarField& Diag = diag();
    Field<Type>& psi =
        const_cast
        <
            GeometricField<Type, faPatchField, areaMesh>&
        >(psi_).primitiveFieldRef();

    forAll(faceLabels, i)
    {
        const label facei = faceLabels[i];
        const Type& value = values[i];

        psi[facei] = value;
        source_[facei] = value*Diag[facei];

        if (symmetric() || asymmetric())
        {
            for (const label edgei : edges[facei])
            {
                if (mesh.isInternalEdge(edgei))
                {
                    if (symmetric())
                    {
                        if (facei == own[edgei])
                        {
                            source_[nei[edgei]] -= upper()[edgei]*value;
                        }
                        else
                        {
                            source_[own[edgei]] -= upper()[edgei]*value;
                        }

                        upper()[edgei] = 0.0;
                    }
                    else
                    {
                        if (facei == own[edgei])
                        {
                            source_[nei[edgei]] -= lower()[edgei]*value;
                        }
                        else
                        {
                            source_[own[edgei]] -= upper()[edgei]*value;
                        }

                        upper()[edgei] = 0.0;
                        lower()[edgei] = 0.0;
                    }
                }
                else
                {
                    const label patchi = mesh.boundary().whichPatch(edgei);

                    if (internalCoeffs_[patchi].size())
                    {
                        const label patchEdgei =
                            mesh.boundary()[patchi].whichEdge(edgei);

                        internalCoeffs_[patchi][patchEdgei] = Zero;
                        boundaryCoeffs_[patchi][patchEdgei] = Zero;
                    }
                }
            }
        }
    }
}


template<class Type>
void Foam::faMatrix<Type>::setValues
(
    const labelUList& faceLabels,
    const Type& value
)
{
    this->setValuesFromList(faceLabels, UniformList<Type>(value));
}


template<class Type>
void Foam::faMatrix<Type>::setValues
(
    const labelUList& faceLabels,
    const UList<Type>& values
)
{
    this->setValuesFromList(faceLabels, values);
}


template<class Type>
void Foam::faMatrix<Type>::setValues
(
    const labelUList& faceLabels,
    const UIndirectList<Type>& values
)
{
    this->setValuesFromList(faceLabels, values);
}


template<class Type>
void Foam::faMatrix<Type>::setReference
(
    const label facei,
    const Type& value,
    const bool forceReference
)
{
    if ((forceReference || psi_.needReference()) && facei >= 0)
    {
        if (Pstream::master())
        {
            source()[facei] += diag()[facei]*value;
            diag()[facei] += diag()[facei];
        }
    }
}


template<class Type>
void Foam::faMatrix<Type>::setReferences
(
    const labelUList& faceLabels,
    const Type& value,
    const bool forceReference
)
{
    if (forceReference || psi_.needReference())
    {
        forAll(faceLabels, facei)
        {
            const label faceId = faceLabels[facei];
            if (faceId >= 0)
            {
                source()[faceId] += diag()[faceId]*value;
                diag()[faceId] += diag()[faceId];
            }
        }
    }
}


template<class Type>
void Foam::faMatrix<Type>::setReferences
(
    const labelUList& faceLabels,
    const UList<Type>& values,
    const bool forceReference
)
{
    if (forceReference || psi_.needReference())
    {
        forAll(faceLabels, facei)
        {
            const label faceId = faceLabels[facei];
            if (faceId >= 0)
            {
                source()[faceId] += diag()[faceId]*values[facei];
                diag()[faceId] += diag()[faceId];
            }
        }
    }
}


template<class Type>
void Foam::faMatrix<Type>::relax(const scalar alpha)
{
    if (alpha <= 0)
    {
        return;
    }

    Field<Type>& S = source();
    scalarField& D = diag();

    // Store the current unrelaxed diagonal for use in updating the source
    scalarField D0(D);

    // Calculate the sum-mag off-diagonal from the interior faces
    scalarField sumOff(D.size(), Zero);
    sumMagOffDiag(sumOff);

    // Handle the boundary contributions to the diagonal
    forAll(psi_.boundaryField(), patchI)
    {
        const faPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (ptf.size())
        {
            const labelUList& pa = lduAddr().patchAddr(patchI);
            Field<Type>& iCoeffs = internalCoeffs_[patchI];

            if (ptf.coupled())
            {
                const Field<Type>& pCoeffs = boundaryCoeffs_[patchI];

                // For coupled boundaries add the diagonal and
                // off-diagonal contributions
                forAll(pa, face)
                {
                    D[pa[face]] += component(iCoeffs[face], 0);
                    sumOff[pa[face]] += mag(component(pCoeffs[face], 0));
                }
            }
            else
            {
                // For non-coupled boundaries subtract the diagonal
                // contribution off-diagonal sum which avoids having to remove
                // it from the diagonal later.
                // Also add the source contribution from the relaxation
                forAll(pa, face)
                {
                    Type iCoeff0 = iCoeffs[face];
                    iCoeffs[face] = cmptMag(iCoeffs[face]);
                    sumOff[pa[face]] -= cmptMin(iCoeffs[face]);
                    iCoeffs[face] /= alpha;
                    S[pa[face]] +=
                        cmptMultiply(iCoeffs[face] - iCoeff0, psi_[pa[face]]);
                }
            }
        }
    }

    // Ensure the matrix is diagonally dominant...
    max(D, D, sumOff);

    // ... then relax
    D /= alpha;

    // Now remove the diagonal contribution from coupled boundaries
    forAll(psi_.boundaryField(), patchI)
    {
        const faPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        if (ptf.size())
        {
            const labelUList& pa = lduAddr().patchAddr(patchI);
            Field<Type>& iCoeffs = internalCoeffs_[patchI];

            if (ptf.coupled())
            {
                forAll(pa, face)
                {
                    D[pa[face]] -= component(iCoeffs[face], 0);
                }
            }
        }
    }

    // Finally add the relaxation contribution to the source.
    S += (D - D0)*psi_.internalField();
}


template<class Type>
void Foam::faMatrix<Type>::relax()
{
    if (psi_.mesh().relaxEquation(psi_.name()))
    {
        relax(psi_.mesh().equationRelaxationFactor(psi_.name()));
    }
    else
    {
        DebugInFunction
            << "Relaxation factor for field " << psi_.name()
            << " not found.  Relaxation will not be used." << endl;
    }
}


template<class Type>
Foam::tmp<Foam::scalarField> Foam::faMatrix<Type>::D() const
{
    tmp<scalarField> tdiag(new scalarField(diag()));
    addCmptAvBoundaryDiag(tdiag.ref());
    return tdiag;
}


template<class Type>
Foam::tmp<Foam::areaScalarField> Foam::faMatrix<Type>::A() const
{
    tmp<areaScalarField> tAphi
    (
        new areaScalarField
        (
            IOobject
            (
                "A("+psi_.name()+')',
                psi_.instance(),
                psi_.db()
            ),
            psi_.mesh(),
            dimensions_/psi_.dimensions()/dimArea,
            zeroGradientFaPatchScalarField::typeName
        )
    );

    tAphi.ref().primitiveFieldRef() = D()/psi_.mesh().S();
    tAphi.ref().correctBoundaryConditions();

    return tAphi;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faPatchField, Foam::areaMesh>>
Foam::faMatrix<Type>::H() const
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> tHphi
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                "H("+psi_.name()+')',
                psi_.instance(),
                psi_.db()
            ),
            psi_.mesh(),
            dimensions_/dimArea,
            zeroGradientFaPatchScalarField::typeName
        )
    );
    GeometricField<Type, faPatchField, areaMesh>& Hphi = tHphi.ref();

    // Loop over field components
    for (direction cmpt=0; cmpt<Type::nComponents; ++cmpt)
    {
        scalarField psiCmpt(psi_.primitiveField().component(cmpt));

        scalarField boundaryDiagCmpt(psi_.size(), Zero);
        addBoundaryDiag(boundaryDiagCmpt, cmpt);
        boundaryDiagCmpt.negate();
        addCmptAvBoundaryDiag(boundaryDiagCmpt);

        Hphi.primitiveFieldRef().replace(cmpt, boundaryDiagCmpt*psiCmpt);
    }

    Hphi.primitiveFieldRef() += lduMatrix::H(psi_.internalField()) + source_;
    addBoundarySource(Hphi.primitiveFieldRef());

    Hphi.primitiveFieldRef() /= psi_.mesh().S();
    Hphi.correctBoundaryConditions();

    return tHphi;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>>
Foam::faMatrix<Type>::flux() const
{
    if (!psi_.mesh().fluxRequired(psi_.name()))
    {
        FatalErrorInFunction
            << "flux requested but " << psi_.name()
            << " not specified in the fluxRequired sub-dictionary of faSchemes"
            << abort(FatalError);
    }

    // construct GeometricField<Type, faePatchField, edgeMesh>
    tmp<GeometricField<Type, faePatchField, edgeMesh>> tfieldFlux
    (
        new GeometricField<Type, faePatchField, edgeMesh>
        (
            IOobject
            (
                "flux("+psi_.name()+')',
                psi_.instance(),
                psi_.db()
            ),
            psi_.mesh(),
            dimensions()
        )
    );
    GeometricField<Type, faePatchField, edgeMesh>& fieldFlux =
        tfieldFlux.ref();

    for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; ++cmpt)
    {
        fieldFlux.primitiveFieldRef().replace
        (
            cmpt,
            lduMatrix::faceH(psi_.primitiveField().component(cmpt))
        );
    }

    FieldField<Field, Type> InternalContrib = internalCoeffs_;

    forAll(InternalContrib, patchI)
    {
        InternalContrib[patchI] =
            cmptMultiply
            (
                InternalContrib[patchI],
                psi_.boundaryField()[patchI].patchInternalField()
            );
    }

    FieldField<Field, Type> NeighbourContrib = boundaryCoeffs_;

    forAll(NeighbourContrib, patchI)
    {
        if (psi_.boundaryField()[patchI].coupled())
        {
            NeighbourContrib[patchI] =
                cmptMultiply
                (
                    NeighbourContrib[patchI],
                    psi_.boundaryField()[patchI].patchNeighbourField()
                );
        }
    }

    forAll(fieldFlux.boundaryField(), patchI)
    {
        fieldFlux.boundaryFieldRef()[patchI] =
            InternalContrib[patchI] - NeighbourContrib[patchI];
    }

    if (faceFluxCorrectionPtr_)
    {
        fieldFlux += *faceFluxCorrectionPtr_;
    }

    return tfieldFlux;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::faMatrix<Type>::operator=(const faMatrix<Type>& famv)
{
    if (this == &famv)
    {
        return;  // Self-assignment is a no-op
    }

    if (&psi_ != &(famv.psi_))
    {
        FatalErrorInFunction
            << "different fields"
            << abort(FatalError);
    }

    lduMatrix::operator=(famv);
    source_ = famv.source_;
    internalCoeffs_ = famv.internalCoeffs_;
    boundaryCoeffs_ = famv.boundaryCoeffs_;

    if (faceFluxCorrectionPtr_ && famv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ = *famv.faceFluxCorrectionPtr_;
    }
    else if (famv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ =
            new GeometricField<Type, faePatchField, edgeMesh>
            (*famv.faceFluxCorrectionPtr_);
    }
}


template<class Type>
void Foam::faMatrix<Type>::operator=(const tmp<faMatrix<Type>>& tfamv)
{
    operator=(tfamv());
    tfamv.clear();
}


template<class Type>
void Foam::faMatrix<Type>::negate()
{
    lduMatrix::negate();
    source_.negate();
    internalCoeffs_.negate();
    boundaryCoeffs_.negate();

    if (faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_->negate();
    }
}


template<class Type>
void Foam::faMatrix<Type>::operator+=(const faMatrix<Type>& famv)
{
    checkMethod(*this, famv, "+=");

    dimensions_ += famv.dimensions_;
    lduMatrix::operator+=(famv);
    source_ += famv.source_;
    internalCoeffs_ += famv.internalCoeffs_;
    boundaryCoeffs_ += famv.boundaryCoeffs_;

    if (faceFluxCorrectionPtr_ && famv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ += *famv.faceFluxCorrectionPtr_;
    }
    else if (famv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ = new
        GeometricField<Type, faePatchField, edgeMesh>
        (
            *famv.faceFluxCorrectionPtr_
        );
    }
}


template<class Type>
void Foam::faMatrix<Type>::operator+=(const tmp<faMatrix<Type>>& tfamv)
{
    operator+=(tfamv());
    tfamv.clear();
}


template<class Type>
void Foam::faMatrix<Type>::operator-=(const faMatrix<Type>& famv)
{
    checkMethod(*this, famv, "+=");

    dimensions_ -= famv.dimensions_;
    lduMatrix::operator-=(famv);
    source_ -= famv.source_;
    internalCoeffs_ -= famv.internalCoeffs_;
    boundaryCoeffs_ -= famv.boundaryCoeffs_;

    if (faceFluxCorrectionPtr_ && famv.faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ -= *famv.faceFluxCorrectionPtr_;
    }
    else if (famv.faceFluxCorrectionPtr_)
    {
        faceFluxCorrectionPtr_ =
            new GeometricField<Type, faePatchField, edgeMesh>
            (-*famv.faceFluxCorrectionPtr_);
    }
}


template<class Type>
void Foam::faMatrix<Type>::operator-=(const tmp<faMatrix<Type>>& tfamv)
{
    operator-=(tfamv());
    tfamv.clear();
}


template<class Type>
void Foam::faMatrix<Type>::operator+=
(
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(*this, su, "+=");
    source() -= su.mesh().S()*su.internalField();
}


template<class Type>
void Foam::faMatrix<Type>::operator+=
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    operator+=(tsu());
    tsu.clear();
}


template<class Type>
void Foam::faMatrix<Type>::operator-=
(
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(*this, su, "-=");
    source() += su.mesh().S()*su.internalField();
}


template<class Type>
void Foam::faMatrix<Type>::operator-=
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    operator-=(tsu());
    tsu.clear();
}


template<class Type>
void Foam::faMatrix<Type>::operator+=
(
    const dimensioned<Type>& su
)
{
    source() -= su.mesh().S()*su;
}


template<class Type>
void Foam::faMatrix<Type>::operator-=
(
    const dimensioned<Type>& su
)
{
    source() += su.mesh().S()*su;
}


template<class Type>
void Foam::faMatrix<Type>::operator*=
(
    const areaScalarField& vsf
)
{
    dimensions_ *= vsf.dimensions();
    lduMatrix::operator*=(vsf.internalField());
    source_ *= vsf.internalField();

    forAll(boundaryCoeffs_, patchI)
    {
        const scalarField psf(vsf.boundaryField()[patchI].patchInternalField());
        internalCoeffs_[patchI] *= psf;
        boundaryCoeffs_[patchI] *= psf;
    }

    if (faceFluxCorrectionPtr_)
    {
        FatalErrorInFunction
            << "cannot scale a matrix containing a faceFluxCorrection"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::faMatrix<Type>::operator*=
(
    const tmp<areaScalarField>& tvsf
)
{
    operator*=(tvsf());
    tvsf.clear();
}


template<class Type>
void Foam::faMatrix<Type>::operator*=
(
    const dimensioned<scalar>& ds
)
{
    dimensions_ *= ds.dimensions();
    lduMatrix::operator*=(ds.value());
    source_ *= ds.value();
    internalCoeffs_ *= ds.value();
    boundaryCoeffs_ *= ds.value();

    if (faceFluxCorrectionPtr_)
    {
        *faceFluxCorrectionPtr_ *= ds.value();
    }
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::checkMethod
(
    const faMatrix<Type>& fam1,
    const faMatrix<Type>& fam2,
    const char* op
)
{
    if (&fam1.psi() != &fam2.psi())
    {
        FatalErrorInFunction
            << "incompatible fields for operation "
            << endl << "    "
            << "[" << fam1.psi().name() << "] "
            << op
            << " [" << fam2.psi().name() << "]"
            << abort(FatalError);
    }

    if
    (
        dimensionSet::checking()
     && fam1.dimensions() != fam2.dimensions()
    )
    {
        FatalErrorInFunction
            << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fam1.psi().name() << fam1.dimensions()/dimArea << " ] "
            << op
            << " [" << fam2.psi().name() << fam2.dimensions()/dimArea << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::checkMethod
(
    const faMatrix<Type>& fam,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const char* op
)
{
    if
    (
        dimensionSet::checking()
     && fam.dimensions()/dimArea != vf.dimensions()
    )
    {
        FatalErrorInFunction
            <<  "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fam.psi().name() << fam.dimensions()/dimArea << " ] "
            << op
            << " [" << vf.name() << vf.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::checkMethod
(
    const faMatrix<Type>& fam,
    const dimensioned<Type>& dt,
    const char* op
)
{
    if
    (
        dimensionSet::checking()
     && fam.dimensions()/dimArea != dt.dimensions()
    )
    {
        FatalErrorInFunction
            << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << fam.psi().name() << fam.dimensions()/dimArea << " ] "
            << op
            << " [" << dt.name() << dt.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
Foam::SolverPerformance<Type> Foam::solve
(
    faMatrix<Type>& fam,
    Istream& solverControls
)
{
    return fam.solve(solverControls);
}


template<class Type>
Foam::SolverPerformance<Type> Foam::solve
(
    const tmp<faMatrix<Type>>& tfam,
    Istream& solverControls
)
{
    SolverPerformance<Type> solverPerf =
        const_cast<faMatrix<Type>&>(tfam()).solve(solverControls);

    tfam.clear();
    return solverPerf;
}


template<class Type>
Foam::SolverPerformance<Type> Foam::solve(faMatrix<Type>& fam)
{
    return fam.solve();
}


template<class Type>
Foam::SolverPerformance<Type> Foam::solve(const tmp<faMatrix<Type>>& tfam)
{
    SolverPerformance<Type> solverPerf =
        const_cast<faMatrix<Type>&>(tfam()).solve();

    tfam.clear();
    return solverPerf;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const faMatrix<Type>& A,
    const faMatrix<Type>& B
)
{
    checkMethod(A, B, "+");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref() += B;
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const tmp<faMatrix<Type>>& tA,
    const faMatrix<Type>& B
)
{
    checkMethod(tA(), B, "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref() += B;
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const faMatrix<Type>& A,
    const tmp<faMatrix<Type>>& tB
)
{
    checkMethod(A, tB(), "+");
    tmp<faMatrix<Type>> tC(tB.ptr());
    tC.ref() += A;
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const tmp<faMatrix<Type>>& tA,
    const tmp<faMatrix<Type>>& tB
)
{
    checkMethod(tA(), tB(), "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref() += tB();
    tB.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const faMatrix<Type>& A
)
{
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().negate();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const tmp<faMatrix<Type>>& tA
)
{
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const faMatrix<Type>& A,
    const faMatrix<Type>& B
)
{
    checkMethod(A, B, "-");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref() -= B;
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const tmp<faMatrix<Type>>& tA,
    const faMatrix<Type>& B
)
{
    checkMethod(tA(), B, "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref() -= B;
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const faMatrix<Type>& A,
    const tmp<faMatrix<Type>>& tB
)
{
    checkMethod(A, tB(), "-");
    tmp<faMatrix<Type>> tC(tB.ptr());
    tC.ref() -= A;
    tC.ref().negate();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const tmp<faMatrix<Type>>& tA,
    const tmp<faMatrix<Type>>& tB
)
{
    checkMethod(tA(), tB(), "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref() -= tB();
    tB.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const faMatrix<Type>& A,
    const faMatrix<Type>& B
)
{
    checkMethod(A, B, "==");
    return (A - B);
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const tmp<faMatrix<Type>>& tA,
    const faMatrix<Type>& B
)
{
    checkMethod(tA(), B, "==");
    return (tA - B);
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const faMatrix<Type>& A,
    const tmp<faMatrix<Type>>& tB
)
{
    checkMethod(A, tB(), "==");
    return (A - tB);
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const tmp<faMatrix<Type>>& tA,
    const tmp<faMatrix<Type>>& tB
)
{
    checkMethod(tA(), tB(), "==");
    return (tA - tB);
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const faMatrix<Type>& A,
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(A, su, "+");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() -= su.mesh().S()*su.internalField();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const tmp<faMatrix<Type>>& tA,
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.mesh().S()*su.internalField();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const faMatrix<Type>& A,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    checkMethod(A, tsu(), "+");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() -= tsu().mesh().S()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const tmp<faMatrix<Type>>& tA,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().S()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const GeometricField<Type, faPatchField, areaMesh>& su,
    const faMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() -= su.mesh().S()*su.internalField();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const GeometricField<Type, faPatchField, areaMesh>& su,
    const tmp<faMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.mesh().S()*su.internalField();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu,
    const faMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "+");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() -= tsu().mesh().S()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu,
    const tmp<faMatrix<Type>>& tA
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= tsu().mesh().S()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const faMatrix<Type>& A,
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(A, su, "-");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() += su.mesh().S()*su.internalField();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const tmp<faMatrix<Type>>& tA,
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += su.mesh().S()*su.internalField();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const faMatrix<Type>& A,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    checkMethod(A, tsu(), "-");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() += tsu().mesh().S()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const tmp<faMatrix<Type>>& tA,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().S()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const GeometricField<Type, faPatchField, areaMesh>& su,
    const faMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().negate();
    tC.ref().source() -= su.mesh().S()*su.internalField();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const GeometricField<Type, faPatchField, areaMesh>& su,
    const tmp<faMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= su.mesh().S()*su.internalField();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu,
    const faMatrix<Type>& A
)
{
    checkMethod(A, tsu(), "-");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().S()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu,
    const tmp<faMatrix<Type>>& tA
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= tsu().mesh().S()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const faMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "+");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() -= su.value()*A.psi().mesh().S();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const tmp<faMatrix<Type>>& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.value()*tC().psi().mesh().S();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const dimensioned<Type>& su,
    const faMatrix<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() -= su.value()*A.psi().mesh().S();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator+
(
    const dimensioned<Type>& su,
    const tmp<faMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() -= su.value()*tC().psi().mesh().S();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const faMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "-");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() += su.value()*tC().psi().mesh().S();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const tmp<faMatrix<Type>>& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += su.value()*tC().psi().mesh().S();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const dimensioned<Type>& su,
    const faMatrix<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().negate();
    tC.ref().source() -= su.value()*A.psi().mesh().S();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator-
(
    const dimensioned<Type>& su,
    const tmp<faMatrix<Type>>& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().negate();
    tC.ref().source() -= su.value()*tC().psi().mesh().S();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const faMatrix<Type>& A,
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(A, su, "==");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() += su.mesh().S()*su.internalField();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const tmp<faMatrix<Type>>& tA,
    const GeometricField<Type, faPatchField, areaMesh>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += su.mesh().S()*su.internalField();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const faMatrix<Type>& A,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    checkMethod(A, tsu(), "==");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() += tsu().mesh().S()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const tmp<faMatrix<Type>>& tA,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tsu
)
{
    checkMethod(tA(), tsu(), "==");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tsu().mesh().S()*tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const faMatrix<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "==");
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref().source() += A.psi().mesh().S()*su.value();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator==
(
    const tmp<faMatrix<Type>>& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref().source() += tC().psi().mesh().S()*su.value();
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator*
(
    const areaScalarField& vsf,
    const faMatrix<Type>& A
)
{
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref() *= vsf;
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator*
(
    const tmp<areaScalarField>& tvsf,
    const faMatrix<Type>& A
)
{
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref() *= tvsf;
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator*
(
    const areaScalarField& vsf,
    const tmp<faMatrix<Type>>& tA
)
{
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref() *= vsf;
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator*
(
    const tmp<areaScalarField>& tvsf,
    const tmp<faMatrix<Type>>& tA
)
{
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref() *= tvsf;
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator*
(
    const dimensioned<scalar>& ds,
    const faMatrix<Type>& A
)
{
    tmp<faMatrix<Type>> tC(new faMatrix<Type>(A));
    tC.ref() *= ds;
    return tC;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::operator*
(
    const dimensioned<scalar>& ds,
    const tmp<faMatrix<Type>>& tA
)
{
    tmp<faMatrix<Type>> tC(tA.ptr());
    tC.ref() *= ds;
    return tC;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const faMatrix<Type>& fam)
{
    os  << static_cast<const lduMatrix&>(fam) << nl
        << fam.dimensions_ << nl
        << fam.source_ << nl
        << fam.internalCoeffs_ << nl
        << fam.boundaryCoeffs_ << endl;

    os.check(FUNCTION_NAME);

    return os;
}


// * * * * * * * * * * * * * * * * Solvers * * * * * * * * * * * * * * * * * //

#include "faMatrixSolve.C"

// ************************************************************************* //
