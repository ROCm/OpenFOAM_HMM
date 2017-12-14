/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "fvMatrix.H"
#include "cellCellStencilObject.H"
#include "oversetFvPatchField.H"

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::dynamicOversetFvMesh::interpolate(Field<T>& psi) const
{
    const cellCellStencil& overlap = Stencil::New(*this);
    const labelListList& stencil = overlap.cellStencil();

    if (stencil.size() != nCells())
    {
        return;
    }

    const mapDistribute& map = overlap.cellInterpolationMap();
    const List<scalarList>& wghts = overlap.cellInterpolationWeights();
    const labelList& cellIDs = overlap.interpolationCells();
    const scalarList& factor = overlap.cellInterpolationWeight();

    Field<T> work(psi);
    map.mapDistributeBase::distribute(work, UPstream::msgType()+1);

    forAll(cellIDs, i)
    {
        label celli = cellIDs[i];

        const scalarList& w = wghts[celli];
        const labelList& nbrs = stencil[celli];
        const scalar f = factor[celli];

        T s(pTraits<T>::zero);
        forAll(nbrs, nbrI)
        {
            s += w[nbrI]*work[nbrs[nbrI]];
        }
        //Pout<< "Interpolated value:" << s << endl;
        //T oldPsi = psi[celli];
        psi[celli] = (1.0-f)*psi[celli] + f*s;
        //Pout<< "psi was:" << oldPsi << " now:" << psi[celli] << endl;
    }
}


template<class GeoField>
void Foam::dynamicOversetFvMesh::interpolate(GeoField& psi) const
{
    interpolate(psi.primitiveFieldRef());
    psi.correctBoundaryConditions();
}


template<class GeoField, class PatchType>
Foam::lduInterfaceFieldPtrsList
Foam::dynamicOversetFvMesh::scalarInterfaces(const GeoField& psi) const
{
    const typename GeoField::Boundary& bpsi = psi.boundaryField();

    lduInterfaceFieldPtrsList interfaces(bpsi.size());

    forAll(interfaces, patchi)
    {
        if
        (
            isA<lduInterfaceField>(bpsi[patchi])
         && isA<PatchType>(bpsi[patchi])
        )
        {
            interfaces.set
            (
                patchi,
                &refCast<const lduInterfaceField>
                (
                    bpsi[patchi]
                )
            );
        }
    }
    return interfaces;
}


template<class Type>
void Foam::dynamicOversetFvMesh::addInterpolation(fvMatrix<Type>& m) const
{
    const cellCellStencilObject& overlap = Stencil::New(*this);
    const List<scalarList>& wghts = overlap.cellInterpolationWeights();
    const labelListList& stencil = overlap.cellStencil();
    const labelList& cellIDs = overlap.interpolationCells();
    const scalarList& factor = overlap.cellInterpolationWeight();
    const labelUList& types = overlap.cellTypes();


    // Force asymmetric matrix (if it wasn't already)
    scalarField& lower = m.lower();
    scalarField& upper = m.upper();
    Field<Type>& source = m.source();
    scalarField& diag = m.diag();


    // Get the addressing. Note that the addressing is now extended with
    // any interpolation faces.
    const lduAddressing& addr = lduAddr();
    const labelUList& upperAddr = addr.upperAddr();
    const labelUList& lowerAddr = addr.lowerAddr();
    const labelUList& ownerStartAddr = addr.ownerStartAddr();
    const labelUList& losortAddr = addr.losortAddr();


    if (!isA<fvMeshPrimitiveLduAddressing>(addr))
    {
        FatalErrorInFunction
            << "Problem : addressing is not fvMeshPrimitiveLduAddressing"
            << exit(FatalError);
    }



    // 1. Adapt lduMatrix for additional faces and new ordering
    upper.setSize(upperAddr.size(), 0.0);
    inplaceReorder(reverseFaceMap_, upper);
    lower.setSize(lowerAddr.size(), 0.0);
    inplaceReorder(reverseFaceMap_, lower);


    // 2. Adapt fvMatrix level: faceFluxCorrectionPtr
    // Question: do we need to do this?
    // This seems to be set/used only by the gaussLaplacianScheme and
    // fvMatrix:correction, both of which are outside the linear solver.



    // Clear out existing connections on cells to be interpolated
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Note: could avoid doing the zeroing of the new faces since these
    //       are set to zero anyway.

    forAll(upperAddr, facei)
    {
        if (types[upperAddr[facei]] == cellCellStencil::INTERPOLATED)
        {
            // Disconnect upper from lower
            label celli = upperAddr[facei];
            lower[facei] *= 1.0-factor[celli];
        }
        if (types[lowerAddr[facei]] == cellCellStencil::INTERPOLATED)
        {
            // Disconnect lower from upper
            label celli = lowerAddr[facei];
            upper[facei] *= 1.0-factor[celli];
        }
    }


    forAll(m.internalCoeffs(), patchI)
    {
        const labelUList& fc = addr.patchAddr(patchI);
        Field<Type>& intCoeffs = m.internalCoeffs()[patchI];
        Field<Type>& bouCoeffs = m.boundaryCoeffs()[patchI];
        forAll(fc, i)
        {
            label celli = fc[i];
            {
                if (types[celli] == cellCellStencil::INTERPOLATED)
                {
                    scalar f = factor[celli];
                    intCoeffs[i] *= 1.0-f;
                    bouCoeffs[i] *= 1.0-f;
                }
                else if (types[celli] == cellCellStencil::HOLE)
                {
                    intCoeffs[i] = pTraits<Type>::zero;
                    bouCoeffs[i] = pTraits<Type>::zero;
                }
            }
        }
    }

    // NOTE: For a non-scalar matrix, in the function fvMatrix::solveSegregated
    // it uses updateInterfaceMatrix to correct for remote source contributions.
    // This unfortunately also triggers the overset interpolation which then
    // produces a non-zero source for interpolated cells.
    // The procedure below calculates this contribution 'correctionSource'
    // and substracts it from the source later in order to compensate.
    Field<Type> correctionSource(diag.size(), pTraits<Type>::zero);

    if (pTraits<Type>::nComponents > 1)
    {
        typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

        typename pTraits<Type>::labelType validComponents
        (
            this->template validComponents<Type>()
        );

        // Get all the overset lduInterfaces
        lduInterfaceFieldPtrsList interfaces
        (
            scalarInterfaces<GeoField, oversetFvPatchField<Type>>
            (
                m.psi()
            )
        );

        const Field<Type>& psiInt = m.psi().primitiveField();

        for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
        {
            if (component(validComponents, cmpt) != -1)
            {
                const scalarField psiCmpt(psiInt.component(cmpt));
                scalarField corrSourceCmpt(correctionSource.component(cmpt));

                forAll(interfaces, patchi)
                {
                    if (interfaces.set(patchi))
                    {
                        interfaces[patchi].initInterfaceMatrixUpdate
                        (
                            corrSourceCmpt,
                            true,
                            psiCmpt,
                            m.boundaryCoeffs()[patchi].component(cmpt),
                            cmpt,
                            Pstream::defaultCommsType
                        );
                    }
                }
                forAll(interfaces, patchi)
                {
                    if (interfaces.set(patchi))
                    {
                        interfaces[patchi].updateInterfaceMatrix
                        (
                            corrSourceCmpt,
                            true,
                            psiCmpt,
                            m.boundaryCoeffs()[patchi].component(cmpt),
                            cmpt,
                            Pstream::defaultCommsType
                        );
                    }
                }
                correctionSource.replace(cmpt, corrSourceCmpt);
            }
        }
    }



    // Modify matrix
    // ~~~~~~~~~~~~~



    // Do hole cells. Note: maybe put into interpolationCells() loop above?
    forAll(types, celli)
    {
        if (types[celli] == cellCellStencil::HOLE)
        {
            //Pout<< "Found hole cell:" << celli
            //    << " at:" << C()[celli]
            //    << " ; setting value" << endl;

            label startLabel = ownerStartAddr[celli];
            label endLabel = ownerStartAddr[celli + 1];

            for (label facei = startLabel; facei < endLabel; facei++)
            {
                upper[facei] = 0.0;
            }

            startLabel = addr.losortStartAddr()[celli];
            endLabel = addr.losortStartAddr()[celli + 1];

            for (label i = startLabel; i < endLabel; i++)
            {
                label facei = losortAddr[i];
                lower[facei] = 0.0;
            }

            const scalar normalisation = V()[celli];
            diag[celli] = normalisation;
            source[celli] = normalisation*m.psi()[celli];
        }
    }


    forAll(cellIDs, i)
    {
        label celli = cellIDs[i];

        const scalar f = factor[celli];
        const scalarList& w = wghts[celli];
        const labelList& nbrs = stencil[celli];
        const labelList& nbrFaces = stencilFaces_[celli];

        const scalar normalisation = V()[celli];

        if (types[celli] == cellCellStencil::HOLE)
        {
            FatalErrorInFunction << "Found HOLE cell " << celli
                << " at:" << C()[celli]
                << " . Should this be in interpolationCells()????"
                << abort(FatalError);
        }
        else
        {
            // Create interpolation stencil. Leave any non-local contributions
            //  to be done by oversetFvPatchField updateMatrixInterface
            // 'correctionSource' is only applied for vector and higher
            // fvMatrix; it is zero for scalar matrices (see above).

            diag[celli] *= (1.0-f);
            diag[celli] += normalisation*f;

            source[celli] *= (1.0-f);
            source[celli] -= correctionSource[celli];

            //Pout<< "Interpolating " << celli << " with factor:" << f
            //    << " new diag:" << diag[celli]
            //    << " new source:" << source[celli]
            //    << " volume:" << V()[celli] << endl;

            forAll(nbrs, nbri)
            {
                label facei = nbrFaces[nbri];

                if (facei != -1)
                {
                    label nbrCelli = nbrs[nbri];

                    // Add the coefficients

                    scalar& u = upper[facei];
                    scalar& l = lower[facei];
                    if (celli < nbrCelli)
                    {
                        u += -normalisation*f*w[nbri];
                    }
                    else
                    {
                        l += -normalisation*f*w[nbri];
                    }
                }
                //else
                //{
                //    Pout<< "addItnerpolation of " << celli
                //        << " at:" << this->cellCentres()[celli]
                //        << " zone:" << this->zoneID()[celli]
                //        << " diag:" << diag[celli]
                //        << " source:" << source[celli]
                //        <<" : skipping remote cell:" << nbrs[nbri]
                //        << " with weight:" << w[nbri]
                //        << endl;
                //}
            }
        }
    }
}


template<class Type>
Foam::SolverPerformance<Type> Foam::dynamicOversetFvMesh::solve
(
    fvMatrix<Type>& m,
    const dictionary& dict
) const
{
    // Check we're running with bcs that can handle implicit matrix manipulation
    const typename GeometricField<Type, fvPatchField, volMesh>::Boundary bpsi =
        m.psi().boundaryField();

    bool hasOverset = false;
    forAll(bpsi, patchi)
    {
        if (isA<oversetFvPatchField<Type>>(bpsi[patchi]))
        {
            hasOverset = true;
            break;
        }
    }

    if (!hasOverset)
    {
        if (debug)
        {
            Pout<< "dynamicOversetFvMesh::solve() :"
                << " bypassing matrix adjustment for field " << m.psi().name()
                << endl;
        }
        return dynamicMotionSolverFvMesh::solve(m, dict);
    }

    if (debug)
    {
        Pout<< "dynamicOversetFvMesh::solve() :"
            << " adjusting matrix for interpolation for field "
            << m.psi().name() << endl;
    }


    // Switch to extended addressing (requires mesh::update() having been
    // called)
    active(true);

    // Adapt matrix
    scalarField oldUpper(m.upper());
    scalarField oldLower(m.lower());
    FieldField<Field, Type> oldInt(m.internalCoeffs());
    FieldField<Field, Type> oldBou(m.boundaryCoeffs());

    addInterpolation(m);

    // Print a bit
    //write(Pout, m, lduAddr());

    // Use lower level solver
    SolverPerformance<Type> s(dynamicMotionSolverFvMesh::solve(m, dict));

    // Restore matrix
    m.upper().transfer(oldUpper);
    m.lower().transfer(oldLower);
    m.internalCoeffs().transfer(oldInt);
    m.boundaryCoeffs().transfer(oldBou);

    // Switch to original addressing
    active(false);

    return s;
}


template<class Type>
void Foam::dynamicOversetFvMesh::write
(
    Ostream& os,
    const fvMatrix<Type>& m,
    const lduAddressing& addr
) const
{
    os  << "** Matrix **" << endl;
    const labelUList& upperAddr = addr.upperAddr();
    const labelUList& lowerAddr = addr.lowerAddr();
    const labelUList& ownerStartAddr = addr.ownerStartAddr();
    const labelUList& losortAddr = addr.losortAddr();

    const scalarField& lower = m.lower();
    const scalarField& upper = m.upper();
    const Field<Type>& source = m.source();
    const scalarField& diag = m.diag();

    forAll(source, celli)
    {
        os  << "cell:" << celli << " diag:" << diag[celli]
            << " source:" << source[celli] << endl;

        label startLabel = ownerStartAddr[celli];
        label endLabel = ownerStartAddr[celli + 1];

        for (label facei = startLabel; facei < endLabel; facei++)
        {
            os  << "    face:" << facei
                << " nbr:" << upperAddr[facei] << " coeff:" << upper[facei]
                << endl;
        }

        startLabel = addr.losortStartAddr()[celli];
        endLabel = addr.losortStartAddr()[celli + 1];

        for (label i = startLabel; i < endLabel; i++)
        {
            label facei = losortAddr[i];
            os  << "    face:" << facei
                << " nbr:" << lowerAddr[facei] << " coeff:" << lower[facei]
                << endl;
        }
    }
    forAll(m.internalCoeffs(), patchi)
    {
        if (m.internalCoeffs().set(patchi))
        {
            const labelUList& fc = addr.patchAddr(patchi);

            os  << "patch:" << patchi << " size:" << fc.size() << endl;
            forAll(fc, i)
            {
                os  << "    cell:" << fc[i]
                    << " int:" << m.internalCoeffs()[patchi][i]
                    << " boun:" << m.boundaryCoeffs()[patchi][i]
                    << endl;
            }
        }
    }


    lduInterfaceFieldPtrsList interfaces =
        m.psi().boundaryField().scalarInterfaces();
    forAll(interfaces, inti)
    {
        if (interfaces.set(inti))
        {
            os  << "interface:" << inti
                << " if type:" << interfaces[inti].interface().type()
                << " fld type:" << interfaces[inti].type() << endl;
        }
    }

    os  << "** End of Matrix **" << endl;
}


template<class GeoField>
void Foam::dynamicOversetFvMesh::correctCoupledBoundaryConditions(GeoField& fld)
{
    typename GeoField::Boundary& bfld = fld.boundaryFieldRef();

    label nReq = Pstream::nRequests();

    forAll(bfld, patchi)
    {
        if (bfld[patchi].coupled())
        {
            Pout<< "initEval of " << bfld[patchi].patch().name() << endl;
            bfld[patchi].initEvaluate(Pstream::defaultCommsType);
        }
    }

    // Block for any outstanding requests
    if (Pstream::parRun())
    {
        Pstream::waitRequests(nReq);
    }

    forAll(bfld, patchi)
    {
        if (bfld[patchi].coupled())
        {
            Pout<< "eval of " << bfld[patchi].patch().name() << endl;
            bfld[patchi].evaluate(Pstream::defaultCommsType);
        }
    }
}


template<class GeoField>
void Foam::dynamicOversetFvMesh::checkCoupledBC(const GeoField& fld)
{
    Pout<< "** starting checking of " << fld.name() << endl;

    GeoField fldCorr(fld.name()+"_correct", fld);
    correctCoupledBoundaryConditions(fldCorr);

    const typename GeoField::Boundary& bfld = fld.boundaryField();
    const typename GeoField::Boundary& bfldCorr = fldCorr.boundaryField();

    forAll(bfld, patchi)
    {
        const typename GeoField::Patch& pfld = bfld[patchi];
        const typename GeoField::Patch& pfldCorr = bfldCorr[patchi];

        Pout<< "Patch:" << pfld.patch().name() << endl;

        forAll(pfld, i)
        {
            if (pfld[i] != pfldCorr[i])
            {
                Pout<< "    " << i << "  orig:" << pfld[i]
                    << " corrected:" << pfldCorr[i] << endl;
            }
        }
    }
    Pout<< "** end of " << fld.name() << endl;
}


// ************************************************************************* //
