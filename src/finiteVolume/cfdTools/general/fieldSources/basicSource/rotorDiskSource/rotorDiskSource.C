/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "rotorDiskSource.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "unitConversion.H"
#include "geometricOneField.H"
#include "fvMatrices.H"
#include "syncTools.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rotorDiskSource, 0);
    addToRunTimeSelectionTable(basicSource, rotorDiskSource, dictionary);

    template<> const char* NamedEnum<rotorDiskSource::geometryModeType, 2>::
        names[] =
    {
        "auto",
        "specified"
    };

    const NamedEnum<rotorDiskSource::geometryModeType, 2>
        rotorDiskSource::geometryModeTypeNames_;

    template<> const char* NamedEnum<rotorDiskSource::inletFlowType, 3>::
        names[] =
    {
        "fixed",
        "surfaceNormal",
        "local"
    };

    const NamedEnum<rotorDiskSource::inletFlowType, 3>
        rotorDiskSource::inletFlowTypeNames_;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::rotorDiskSource::checkData()
{
    switch (selectionMode())
    {
        case smCellSet:
        case smCellZone:
        case smAll:
        {
            // set the profile ID for each blade section
            profiles_.connectBlades(blade_.profileName(), blade_.profileID());
            switch (inletFlow_)
            {
                case ifFixed:
                {
                    coeffs_.lookup("inletVelocity") >> inletVelocity_;
                    break;
                }
                case ifSurfaceNormal:
                {
                    scalar UIn
                    (
                        readScalar(coeffs_.lookup("inletNormalVelocity"))
                    );
                    inletVelocity_ = -coordSys_.e3()*UIn;
                    break;
                }
                case ifLocal:
                {
                    // do nothing
                    break;
                }
                default:
                {
                    FatalErrorIn("void rotorDiskSource::checkData()")
                        << "Unknown inlet velocity type" << abort(FatalError);
                }
            }


            break;
        }
        default:
        {
            FatalErrorIn("void rotorDiskSource::checkData()")
                << "Source cannot be used with '"
                << selectionModeTypeNames_[selectionMode()]
                << "' mode.  Please use one of: " << nl
                << selectionModeTypeNames_[smCellSet] << nl
                << selectionModeTypeNames_[smCellZone] << nl
                << selectionModeTypeNames_[smAll]
                << exit(FatalError);
        }
    }
}


void Foam::rotorDiskSource::setFaceArea(vector& axis, const bool correct)
{
    static const scalar tol = 0.8;

    const label nInternalFaces = mesh_.nInternalFaces();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const vectorField& Sf = mesh_.Sf().internalField();
    const scalarField& magSf = mesh_.magSf().internalField();

    vector n = vector::zero;
    label nFace = 0;

    // calculate cell addressing for selected cells
    labelList cellAddr(mesh_.nCells(), -1);
    UIndirectList<label>(cellAddr, cells_) = identity(cells_.size());
    labelList nbrFaceCellAddr(mesh_.nFaces() - nInternalFaces, -1);
    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label faceI = pp.start() + i;
                label nbrFaceI = faceI - nInternalFaces;
                label own = mesh_.faceOwner()[faceI];
                nbrFaceCellAddr[nbrFaceI] = cellAddr[own];
            }
        }
    }

    // correct for parallel running
    syncTools::swapBoundaryFaceList(mesh_, nbrFaceCellAddr);


    // add internal field contributions
    for (label faceI = 0; faceI < nInternalFaces; faceI++)
    {
        const label own = cellAddr[mesh_.faceOwner()[faceI]];
        const label nbr = cellAddr[mesh_.faceNeighbour()[faceI]];

        if ((own != -1) && (nbr == -1))
        {
            vector nf = Sf[faceI]/magSf[faceI];

            if ((nf & axis) > tol)
            {
                area_[own] += magSf[faceI];
                n += Sf[faceI];
                nFace++;
            }
        }
        else if ((own == -1) && (nbr != -1))
        {
            vector nf = Sf[faceI]/magSf[faceI];

            if ((-nf & axis) > tol)
            {
                area_[nbr] += magSf[faceI];
                n -= Sf[faceI];
                nFace++;
            }
        }
    }


    // add boundary contributions
    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];
        const vectorField& Sfp = mesh_.Sf().boundaryField()[patchI];
        const scalarField& magSfp = mesh_.magSf().boundaryField()[patchI];

        if (pp.coupled())
        {
            forAll(pp, j)
            {
                const label faceI = pp.start() + j;
                const label own = cellAddr[mesh_.faceOwner()[faceI]];
                const bool nbr = nbrFaceCellAddr[faceI - nInternalFaces];
                const vector nf = Sfp[j]/magSfp[j];

                if ((own != -1) && (nbr == -1) && ((nf & axis) > tol))
                {
                    area_[own] += magSfp[j];
                    n += Sfp[j];
                    nFace++;
                }
            }
        }
        else
        {
            forAll(pp, j)
            {
                const label faceI = pp.start() + j;
                const label own = cellAddr[mesh_.faceOwner()[faceI]];
                const vector nf = Sfp[j]/magSfp[j];

                if ((own != -1) && ((nf & axis) > 0.8))
                {
                    area_[own] += magSfp[j];
                    n += Sfp[j];
                    nFace++;
                }
            }
        }
    }

    if (correct)
    {
        reduce(n, sumOp<vector>());
        axis = n/mag(n);
    }
}


void Foam::rotorDiskSource::createCoordinateSystem()
{
    // construct the local rotor co-prdinate system
    vector origin(vector::zero);
    vector axis(vector::zero);
    vector refDir(vector::zero);

    geometryModeType gm =
        geometryModeTypeNames_.read(coeffs_.lookup("geometryMode"));

    switch (gm)
    {
        case gmAuto:
        {
            // determine rotation origin
            scalar sumV = 0.0;
            const scalarField& V = mesh_.V();
            const vectorField& C = mesh_.C();
            forAll(cells_, i)
            {
                const label cellI = cells_[i];
                sumV += V[cellI];
                origin += V[cellI]*C[cellI];
            }
            origin /= sumV;

            // determine first radial vector
            vector dx1(vector::zero);
            scalar magR = -GREAT;
            forAll(cells_, i)
            {
                const label cellI = cells_[i];
                vector test = C[cellI] - origin;
                if (mag(test) > magR)
                {
                    dx1 = test;
                    magR = mag(test);
                }
            }

            // determine second radial vector and cross to determine axis
            forAll(cells_, i)
            {
                const label cellI = cells_[i];
                vector dx2 = C[cellI] - origin;
                if (mag(dx2) > 0.5*magR)
                {
                    axis = dx1 ^ dx2;
                    if (mag(axis) > SMALL)
                    {
                        break;
                    }
                }
            }
            axis /= mag(axis);

            // axis direction is somewhat arbitrary - check if user needs
            // needs to reverse
            bool reverse(readBool(coeffs_.lookup("reverseAxis")));
            if (reverse)
            {
                axis *= -1.0;
            }

            coeffs_.lookup("refDirection") >> refDir;

            // set the face areas and apply correction to calculated axis
            // e.g. if cellZone is more than a single layer in thickness
            setFaceArea(axis, true);

            break;
        }
        case gmSpecified:
        {
            coeffs_.lookup("origin") >> origin;
            coeffs_.lookup("axis") >> axis;
            coeffs_.lookup("refDirection") >> refDir;

            setFaceArea(axis, false);

            break;
        }
        default:
        {
            FatalErrorIn("rotorDiskSource::createCoordinateSystem()")
                << "Unknown geometryMode " << geometryModeTypeNames_[gm]
                << ". Available geometry modes include "
                << geometryModeTypeNames_ << exit(FatalError);
        }
    }

    coordSys_ = cylindricalCS("rotorCoordSys", origin, axis, refDir, false);

    const scalar sumArea = gSum(area_);
    const scalar diameter = Foam::sqrt(4.0*sumArea/mathematical::pi);
    Info<< "    Rotor gometry:" << nl
        << "    - disk diameter = " << diameter << nl
        << "    - disk area     = " << sumArea << nl
        << "    - origin        = " << coordSys_.origin() << nl
        << "    - r-axis        = " << coordSys_.e1() << nl
        << "    - psi-axis      = " << coordSys_.e2() << nl
        << "    - z-axis        = " << coordSys_.e3() << endl;
}


void Foam::rotorDiskSource::constructGeometry()
{
    const vectorField& C = mesh_.C();

    const vector rDir = coordSys_.e1();
    const vector zDir = coordSys_.e3();

    forAll(cells_, i)
    {
        const label cellI = cells_[i];

        // position in rotor co-ordinate system
        x_[i] = coordSys_.localPosition(C[cellI]);

        // cache max radius
        rMax_ = max(rMax_, x_[i].x());

        // determine swept angle relative to rDir axis
        scalar psi = x_[i].y() - rDir.y();

        // blade flap angle
        scalar beta = flap_.beta0 - flap_.beta1*cos(psi) - flap_.beta2*sin(psi);

        // determine rotation tensor to convert into the rotor cone plane
        scalar c = cos(-beta);
        scalar s = sin(-beta);
        R_[i] = tensor(1, 0, 0, 0, c, s, 0, -s, c);

        // geometric angle of attack - not including twist
        alphag_[i] = trim_.alphaC - trim_.A*cos(psi) - trim_.B*sin(psi);
    }
}


Foam::tmp<Foam::vectorField> Foam::rotorDiskSource::inflowVelocity
(
    const volVectorField& U
) const
{
    switch (inletFlow_)
    {
        case ifFixed:
        case ifSurfaceNormal:
        {
            return tmp<vectorField>
            (
                new vectorField(mesh_.nCells(), inletVelocity_)
            );

            break;
        }
        case ifLocal:
        {
            return U.internalField();

            break;
        }
        default:
        {
            FatalErrorIn
            (
                "Foam::tmp<Foam::vectorField> "
                "Foam::rotorDiskSource::inflowVelocity"
                "(const volVectorField&) const"
            )   << "Unknown inlet flow specification" << abort(FatalError);
        }
    }

    return tmp<vectorField>(new vectorField(mesh_.nCells(), vector::zero));
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::rotorDiskSource::rotorDiskSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh

)
:
    basicSource(name, modelType, dict, mesh),
    rhoName_("none"),
    omega_(0.0),
    nBlades_(0),
    inletFlow_(ifLocal),
    inletVelocity_(vector::zero),
    tipEffect_(1.0),
    flap_(),
    trim_(),
    blade_(coeffs_.subDict("blade")),
    profiles_(coeffs_.subDict("profiles")),
    x_(cells_.size(), vector::zero),
    R_(cells_.size(), I),
    alphag_(cells_.size(), 0.0),
    area_(cells_.size(), 0.0),
    coordSys_(false),
    rMax_(0.0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rotorDiskSource::~rotorDiskSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rotorDiskSource::addSup(fvMatrix<vector>& eqn, const label fieldI)
{
    // add source to rhs of eqn

    const volVectorField& U = eqn.psi();

    if (eqn.dimensions() == dimForce)
    {
        coeffs_.lookup("rhoName") >> rhoName_;

        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        eqn -= calculateForces
            (
                rho.internalField(),
                inflowVelocity(U),
                dimForce/dimVolume
            );
    }
    else
    {
        eqn -= calculateForces
            (
                oneField(),
                inflowVelocity(U),
                dimForce/dimVolume/dimDensity
            );
    }
}


void Foam::rotorDiskSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::rotorDiskSource::read(const dictionary& dict)
{
    if (basicSource::read(dict))
    {
        coeffs_.lookup("fieldNames") >> fieldNames_;
        applied_.setSize(fieldNames_.size(), false);

        scalar rpm(readScalar(coeffs_.lookup("rpm")));
        omega_ = rpm/60.0*mathematical::twoPi;

        coeffs_.lookup("nBlades") >> nBlades_;

        inletFlow_ = inletFlowTypeNames_.read(coeffs_.lookup("inletFlowType"));

        coeffs_.lookup("tipEffect") >> tipEffect_;

        const dictionary& flapCoeffs(coeffs_.subDict("flapCoeffs"));
        flapCoeffs.lookup("beta0") >> flap_.beta0;
        flapCoeffs.lookup("beta1") >> flap_.beta1;
        flapCoeffs.lookup("beta2") >> flap_.beta2;
        flap_.beta0 = degToRad(flap_.beta0);

        const dictionary& trimCoeffs(coeffs_.subDict("trimCoeffs"));
        trimCoeffs.lookup("alphaC") >> trim_.alphaC;
        trimCoeffs.lookup("A") >> trim_.A;
        trimCoeffs.lookup("B") >> trim_.B;
        trim_.alphaC = degToRad(trim_.alphaC);

        checkData();

        createCoordinateSystem();

        constructGeometry();

        if (debug)
        {
            writeField("alphag", alphag_, true);
            writeField("faceArea", area_, true);
        }

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
