/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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

#include "updateMethod.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(updateMethod, 0);
    defineRunTimeSelectionTable(updateMethod, dictionary);
}


// * * * * * * * * * *  Protected  Member Functions  * * * * * * * * * * * * //

const Foam::scalarField Foam::updateMethod::leftMult
(
    const scalarField& s,
    const SquareMatrix<scalar>& m
)
{
    if (s.size() != m.n())
    {
        FatalErrorInFunction
            << "scalar derivative and HessianInv matrix do not have the "
            << "same dimension"
            << abort(FatalError);
    }

    scalarField res(s.size(), Zero);
    forAll(s, i)
    {
        forAll(s, j)
        {
            res[i] += s[j]*m[j][i];
        }
    }

    return (res);
}


const Foam::scalarField Foam::updateMethod::rightMult
(
    const SquareMatrix<scalar>& m,
    const scalarField& s
)
{
    if (s.size() != m.n())
    {
        FatalErrorInFunction
            << "scalar derivative and HessianInv matrix do not have the "
            << "same dimension"
            << abort(FatalError);
    }

    scalarField res(s.size(), Zero);
    forAll(s, i)
    {
        forAll(s, j)
        {
            res[i] += m[i][j]*s[j];
        }
    }

    return (res);
}


Foam::SquareMatrix<Foam::scalar> Foam::updateMethod::outerProd
(
    const scalarField& a,
    const scalarField& b
)
{
    if (a.size() != b.size())
    {
        FatalErrorInFunction
            << "operands of outerProduct do not have the same dimension"
            << abort(FatalError);
    }

    SquareMatrix<scalar> res(a.size(), Zero);
    forAll(a, i)
    {
        forAll(a, j)
        {
            res[i][j]  = a[i]*b[j];
        }
    }

    return (res);
}


Foam::SquareMatrix<Foam::scalar>
Foam::updateMethod::inv(SquareMatrix<scalar> A)
{
    label n(A.n());
    SquareMatrix<scalar> invA(n, Zero);

    // LU decomposition of A
    labelList pivotIndices(n, Zero);
    LUDecompose(A, pivotIndices);
    DebugInfo
        << "LU decomposed A " << A << endl;

    // Compute inverse of A by successive back-substitutions.
    for (label j = 0; j < n; j++)
    {
        scalarField rhs(n, Zero);
        rhs[j] = scalar(1);
        LUBacksubstitute(A, pivotIndices, rhs);
        // After LUBacksubstitute, rhs contains the j-th column of the inverse
        for (label i = 0; i < n; i++)
        {
            invA[i][j] = rhs[i];
        }
    }


    /*
    // Alternative using SVD. Slower and less accurate
    tempscalarRectangularMatrix Atemp(n, n, 0);
    for (label i = 0; i < n; i++)
    {
        for (label j = 0; j < n; j++)
        {
            Atemp[i][j] = A[i][j];
        }
    }
    scalarRectangularMatrix invTemp = SVDinv(Atemp);
    scalarSquareMatrix invA(n, n, 0);
    for (label i = 0; i < n; i++)
    {
        for (label j = 0; j < n; j++)
        {
            invA[i][j] = invTemp[i][j];
        }
    }
    */

    return invA;
}


Foam::scalar Foam::updateMethod::globalSum(const scalarField& field)
{
    scalar value(0);
    if (globalSum_)
    {
        value = gSum(field);
    }
    else
    {
        value = sum(field);
    }
    return value;
}


Foam::scalar Foam::updateMethod::globalSum(tmp<scalarField>& tfield)
{
    scalar value = globalSum(tfield());
    tfield.clear();
    return value;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::updateMethod::updateMethod
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    optMethodIODict_
    (
        IOobject
        (
            "updateMethodDict",
            mesh_.time().timeName(),
            "uniform",
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    objectiveDerivatives_(0),
    constraintDerivatives_(0),
    objectiveValue_(0),
    cValues_(0),
    correction_(0),
    cumulativeCorrection_(0),
    eta_(1),
    initialEtaSet_(false),
    correctionFolder_(mesh_.time().globalPath()/"optimisation"/"correction"),
    globalSum_
    (
        dict.getOrDefault<bool>("globalSum", false)
    )
{
    // Create folder to store corrections
    if (Pstream::master())
    {
        mkDir(correctionFolder_);
    }

    // Set initial eta, if present. It might be set either in the
    // optimisationDict or in the specific dictionary dedicated to the
    // updateMethod
    if (dict.readIfPresent("eta", eta_))
    {
        initialEtaSet_ = true;
    }
    else if (optMethodIODict_.readIfPresent("eta", eta_))
    {
        initialEtaSet_ = true;
    }
}


Foam::dictionary Foam::updateMethod::coeffsDict()
{
    return dict_.subOrEmptyDict(type());
}


// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::updateMethod> Foam::updateMethod::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    const word modelType(dict.get<word>("method"));

    Info<< "updateMethod type : " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "updateMethod",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<updateMethod>(ctorPtr(mesh, dict));
}


// * * * * * * * * * * * * * * *  Member Functions   * * * * * * * * * * * * //

void Foam::updateMethod::setObjectiveDeriv(const scalarField& derivs)
{
    objectiveDerivatives_ = derivs;
}


void Foam::updateMethod::setConstraintDeriv
(
    const PtrList<scalarField>& derivs
)
{
    constraintDerivatives_ = derivs;
}


void Foam::updateMethod::setObjectiveValue(const scalar value)
{
    objectiveValue_ = value;
}


void Foam::updateMethod::setConstraintValues(const scalarField& values)
{
    cValues_ = values;
}


void Foam::updateMethod::setStep(const scalar eta)
{
    eta_ = eta;
}


void Foam::updateMethod::setGlobalSum(const bool useGlobalSum)
{
    globalSum_ = useGlobalSum;
}


Foam::scalarField& Foam::updateMethod::returnCorrection()
{
    computeCorrection();
    return correction_;
}


void Foam::updateMethod::writeCorrection()
{
    if (Pstream::master())
    {
        // Allocate cumulativeCorrection if necessary
        if (cumulativeCorrection_.empty())
        {
            cumulativeCorrection_.setSize(correction_.size(), Zero);
        }
        // Accumulate correction
        cumulativeCorrection_ += correction_;

        fileName correctionFile
        (
            correctionFolder_/"correction" + mesh_.time().timeName()
        );
        fileName cumulativeCorrectionFile
        (
            correctionFolder_/"cumulativeCorrection" + mesh_.time().timeName()
        );

        OFstream corFile(correctionFile);
        OFstream cumulCorFile(cumulativeCorrectionFile);
        forAll(correction_, cI)
        {
            corFile
                << cI << " " << correction_[cI] << endl;
            cumulCorFile
                << cI << " " << cumulativeCorrection_[cI] << endl;
        }
    }
}


Foam::scalar Foam::updateMethod::computeMeritFunction()
{
    return objectiveValue_;
}


Foam::scalar Foam::updateMethod::meritFunctionDirectionalDerivative()
{
    return globalSum(objectiveDerivatives_*correction_);
}


bool& Foam::updateMethod::initialEtaSet()
{
    return initialEtaSet_;
}


void Foam::updateMethod::updateOldCorrection
(
    const scalarField& oldCorrection
)
{
    correction_ = oldCorrection;
}


void Foam::updateMethod::write()
{
    // Insert eta if set
    if (initialEtaSet_)
    {
        optMethodIODict_.add<scalar>("eta", eta_, true);
    }

    optMethodIODict_.add<scalarField>("correction", correction_, true);

    // Write IOdictionary
    // Always write in ASCII format.
    // Even when choosing to write in binary through controlDict,
    // the content is written in ASCII format but with a binary header.
    // This creates problems when the content is read back in
    // (e.g. continuation)
    optMethodIODict_.regIOobject::writeObject
    (
        IOstreamOption(IOstream::ASCII, mesh_.time().writeCompression()),
        true
    );
}


// ************************************************************************* //
