/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd
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

#include "laserDTRM.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "interpolationCell.H"
#include "interpolationCellPoint.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(laserDTRM, 0);
        addToRadiationRunTimeSelectionTables(laserDTRM);
    }

    defineTemplateTypeNameAndDebugWithName
    (
        Cloud<DTRMParticle>,
        "DTRMCloud",
        0
    );

    namespace radiation
    {
        template<>
        const char* Foam::NamedEnum
        <
            Foam::radiation::laserDTRM::powerDistributionMode,
            3
        >::names[] =
        {
            "Gaussian",
            "manual",
            "uniform"
        };
    }
}

const Foam::NamedEnum
<
    Foam::radiation::laserDTRM::powerDistributionMode,
    3
> Foam::radiation::laserDTRM::powerDistypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::radiation::laserDTRM::calculateIp(scalar r, scalar theta)
{
    const scalar t = mesh_.time().value();
    switch(mode_)
    {
        case pdGaussian:
        {
            scalar I0 =
                laserPower_->value(t)/(mathematical::twoPi*sqr(sigma_));

            return(I0*exp(-sqr(r)/2.0/sqr(sigma_)));

            break;
        }
        case pdManual:
        {
            return(laserPower_->value(t)*powerDistribution_()(theta, r));
            break;
        }
        case pdUniform:
        {
            return
            (
                laserPower_->value(t)/(mathematical::pi*sqr(focalLaserRadius_))
            );
        }
        default:
        {
            FatalErrorInFunction
                    << "Unhandled type " << powerDistypeNames_
                    << abort(FatalError);
            return(0);
        }
    }
}


Foam::tmp<Foam::volVectorField> Foam::radiation::laserDTRM::nHatfv
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    const dimensionedScalar deltaN
    (
        "deltaN",
        1e-8/pow(average(mesh_.V()), 1.0/3.0)
    );

    const dimensionedScalar minAlpha
    (
        "minAlpha", dimless, 1e-3
    );

    volVectorField gradAlphaf
    (
        "gradAlphaf",
         (alpha2 + minAlpha)*fvc::grad(alpha1)
       - (alpha1 + minAlpha)*fvc::grad(alpha2)
    );

   // Face unit interface normal
   return gradAlphaf/(mag(gradAlphaf)+ deltaN);
}


Foam::tmp<Foam::volScalarField>Foam::radiation::laserDTRM::nearInterface
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    return
        pos(alpha1 - 0.1)*pos(0.9 - alpha1)
      * pos(alpha2 - 0.1)*pos(0.9 - alpha2);
}


void Foam::radiation::laserDTRM::initialise()
{
    // Initialise the DTRM particles
    DTRMCloud_.clear();

    const scalar t = mesh_.time().value();
    const vector lPosition = focalLaserPosition_->value(t);
    vector lDir = laserDirection_->value(t);
    lDir /= mag(lDir);

    if (debug)
    {
        Info << "Laser position : " << lPosition << endl;
        Info << "Laser direction : " << lDir << endl;
    }

    // Find a vector on the area plane. Normal to laser direction
    vector rArea = vector::zero;
    scalar magr = 0.0;
    cachedRandom rnd(1234, -1);
    while (magr < VSMALL)
    {
        vector v = rnd.sample01<vector>();
        rArea = v - (v & lDir)*lDir;
        magr = mag(rArea);
    }
    rArea /= mag(rArea);

    scalar dr =  focalLaserRadius_/ndr_;
    scalar dTheta =  mathematical::twoPi/ndTheta_;

    nParticles_ = ndr_*ndTheta_;

    switch(mode_)
    {
        case pdGaussian:
        {
            sigma_ = readScalar(lookup("sigma"));
            break;
        }
        case pdManual:
        {
            powerDistribution_.reset
            (
                new interpolation2DTable<scalar>(*this)
            );

            // Check dimensions ndr and ndTheta
//             if
//             (
//                (powerDistribution_->size() != ndTheta_)
//             || (powerDistribution_().first().second().size() != ndr_)
//             )
//             {
//                  FatalErrorInFunction
//                     << " The table dimensions should correspond with ndTheta "
//                     << " and ndr "
//                     << exit(FatalError);
//             }

            break;
        }
        case pdUniform:
        {
            break;
        }
    }

    // Target position
    point p1 = vector::zero;

    // Seed DTRM particles
    // TODO: currently only applicable to 3-D cases
    point p0 = lPosition;
    scalar power(0.0);
    scalar area(0.0);
    p1 = p0;
    if (mesh_.nGeometricD() == 3)
    {
        //scalar r0 = dr/2.0;
        //scalar r1Max0 = drMax/2.0;

        for (label ri = 0; ri < ndr_; ri++)
        {
            scalar r1 = SMALL + dr*ri;

            scalar r2 = r1 + dr;

            scalar rP = ((r1 + r2)/2);

            // local radius on disk
            vector localR = ((r1 + r2)/2)*rArea;

            // local final radius on disk
            vector finalR = rP*rArea;

            scalar theta0 = 0.0;//dTheta/2.0;
            for (label thetai = 0; thetai < ndTheta_; thetai++)
            {
                scalar theta1 = theta0 + SMALL  + dTheta*thetai;

                scalar theta2 = theta1 + dTheta;

                scalar thetaP = (theta1 + theta2)/2.0;

                quaternion Q(lDir, thetaP);

                // Initial position on disk
                vector initialPos = (Q.R() & localR);

                // Final position on disk
                vector finalPos = (Q.R() & finalR);

                // Initial position
                p0 = lPosition + initialPos;

                // calculate target point using new deviation rl
                p1 = lPosition + finalPos + (0.5*maxTrackLength_*lDir);

                //scalar p = magSqr(p0 - lPosition);

                scalar Ip = calculateIp(rP, thetaP);

                scalar dAi = (sqr(r2) - sqr(r1))*(theta2 - theta1)/2.0;

                power += Ip*dAi;
                area += dAi;

                label cellI = mesh_.findCell(p0);

                if (cellI != -1)
                {
                    // Create a new particle
                    DTRMParticle* pPtr = new DTRMParticle
                        (mesh_, p0, p1, Ip, cellI, dAi, 0 , 0.01*Ip, true);

                    // Add to cloud
                    DTRMCloud_.addParticle(pPtr);
                }

                if (returnReduce(cellI, maxOp<label>()) == -1)
                {
                    WarningIn("void Foam::radiation::laserDTRM::initialise()")
                        << "Cannot find owner cell for particle at position " << p0
                        << endl;
                }
            }
        }
    }
    else
    {
            FatalErrorIn("void Foam::radiation::laserDTRM::initialise()")
            << "Current functionality limited to 3-D cases"
            << exit(FatalError);
    }

    if (debug)
    {
        Info << "Total Power in the laser : " << power << endl;
        Info << "Total Area in the laser : " << area << endl;
        Info << "Number of particles in the laser : "
             << returnReduce(DTRMCloud_.size(), sumOp<label>()) << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::laserDTRM::laserDTRM(const volScalarField& T)
:
    radiationModel(typeName, T),
    mode_(powerDistypeNames_.read(lookup("mode"))),
    DTRMCloud_(mesh_, "DTRMCloud", IDLList<DTRMParticle>()),
    nParticles_(0),
    ndTheta_(readLabel(lookup("nTheta"))),
    ndr_(readLabel(lookup("nr"))),
    maxTrackLength_(mesh_.bounds().mag()),

    focalLaserPosition_
    (
        Function1<point>::New("focalLaserPosition", *this)
    ),
    laserDirection_
    (
        Function1<vector>::New("laserDirection", *this)
    ),


    focalLaserRadius_(readScalar(lookup("focalLaserRadius"))),
    qualityBeamLaser_
    (
        lookupOrDefault<scalar>("qualityBeamLaser", 0.0)
    ),

    sigma_(0),
    laserPower_(Function1<scalar>::New("laserPower", *this)),
    powerDistribution_(),

    reflection_(reflectionModel::New(*this, mesh_).ptr()),
    reflectionSwitch_(lookupOrDefault("reflection", false)),
    initialPhase_(lookupOrDefault("initialPhase", word::null)),
    alpha1_(lookupOrDefault("alpha1", word::null)),
    alpha2_(lookupOrDefault("alpha2", word::null)),
    Qin_
    (
        IOobject
        (
            "Qin",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qin", dimPower/dimArea, 0.0)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
    ),
    Q_
    (
     IOobject
        (
            "Q",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Q", dimPower/dimVolume, 0.0)
    )
{
    initialise();
}


Foam::radiation::laserDTRM::laserDTRM
(
    const dictionary& dict,
    const volScalarField& T
)
:
    radiationModel(typeName, dict, T),
    mode_(powerDistypeNames_.read(lookup("mode"))),
    DTRMCloud_(mesh_, "DTRMCloud", IDLList<DTRMParticle>()),
    nParticles_(0),
    ndTheta_(readLabel(lookup("nTheta"))),
    ndr_(readLabel(lookup("nr"))),
    maxTrackLength_(mesh_.bounds().mag()),

    focalLaserPosition_
    (
        Function1<point>::New("focalLaserPosition", *this)
    ),
    laserDirection_
    (
        Function1<vector>::New("laserDirection", *this)
    ),

    focalLaserRadius_(readScalar(lookup("focalLaserRadius"))),
    qualityBeamLaser_
    (
        lookupOrDefault<scalar>("qualityBeamLaser", 0.0)
    ),

    sigma_(0),
    laserPower_(Function1<scalar>::New("laserPower", *this)),
    powerDistribution_(),

    reflection_(reflectionModel::New(*this, mesh_).ptr()),
    reflectionSwitch_(dict.lookupOrDefault("reflection", false)),
    initialPhase_(lookupOrDefault("initialPhase", word::null)),
    alpha1_(lookupOrDefault("alpha1", word::null)),
    alpha2_(lookupOrDefault("alpha2", word::null)),
    Qin_
    (
        IOobject
        (
            "Qin",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qin", dimPower/dimArea, 0.0)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
    ),
    Q_
    (
     IOobject
        (
            "Q",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Q", dimPower/pow3(dimLength), 0.0)
    )
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::laserDTRM::~laserDTRM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::laserDTRM::read()
{
    if (radiationModel::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiation::laserDTRM::calculate()
{
    tmp<volScalarField> treflectingCells
    (
        new volScalarField
        (
            IOobject
            (
                "reflectingCellsVol",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0)
        )
    );
    volScalarField& reflectingCellsVol = treflectingCells.ref();

    // Reset the fields
    Qin_ == dimensionedScalar("zero", Qin_.dimensions(), 0);
    Q_ == dimensionedScalar("zero", Q_.dimensions(), 0);

    a_ = absorptionEmission_->a();
    e_ = absorptionEmission_->e();
    E_ = absorptionEmission_->E();

    const interpolationCell<scalar> aInterp(a_);
    const interpolationCell<scalar> eInterp(e_);
    const interpolationCell<scalar> EInterp(E_);
    const interpolationCell<scalar> TInterp(T_);

    labelField reflectingCells(mesh_.nCells(), 0);

    tmp<volVectorField> tnHat
    (
        new volVectorField
        (
            IOobject
            (
                "nHat",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimless, vector::zero)
        )
    );
    volVectorField& nHat = tnHat.ref();

    autoPtr<interpolationCellPoint<vector>> nHatIntrPtr;

    if (reflectionSwitch_)
    {
        const volScalarField& initialPhase =
            mesh_.lookupObject<volScalarField>(initialPhase_);

        if (alpha1_ != word::null)
        {
            const volScalarField& alpha1 =
                mesh_.lookupObject<volScalarField>(alpha1_);

            nHat = nHatfv(initialPhase, alpha1);

            forAll(alpha1, cellI)
            {
                if ((alpha1[cellI] > 0.9) && (mag(nHat[cellI]) > 0.9))
                {
                    reflectingCells[cellI] = 1;
                    reflectingCellsVol[cellI] = 1.0;
                }
            }
        }

        if (alpha2_ != word::null)
        {
            const volScalarField& alpha2 =
                mesh_.lookupObject<volScalarField>(alpha2_);

            nHat += nHatfv(initialPhase, alpha2);

            forAll(alpha2, cellI)
            {
                if ((alpha2[cellI] > 0.9) && (mag(nHat[cellI]) > 0.9))
                {
                    reflectingCells[cellI] = 1;
                    reflectingCellsVol[cellI] = 1.0;
                }
            }
        }

        nHatIntrPtr.reset
        (
            new interpolationCellPoint<vector>(nHat)
        );

    }

    DTRMParticle::trackingData td
    (
        DTRMCloud_,
        aInterp,
        eInterp,
        EInterp,
        TInterp,
        nHatIntrPtr,
        reflectingCells,
        reflection_,
        Q_
    );

    Info << "Move particles..."
         << returnReduce(DTRMCloud_.size(), sumOp<label>()) << endl;

    DTRMCloud_.move(td, mesh_.time().deltaTValue());

    // Normalize by cell volume
    Q_.primitiveFieldRef() /= mesh_.V();

    if (debug)
    {

        OFstream osRef
        (
            type() + ":particlePath.obj"
        );
        label vertI = 0;

        forAllIter(Cloud<DTRMParticle>, DTRMCloud_, iter)
        {
            DTRMParticle& p = iter();
            meshTools::writeOBJ(osRef, p.position());
            vertI++;
            meshTools::writeOBJ(osRef, p.p0());
            vertI++;
            osRef << "l " << vertI-1 << ' ' << vertI << nl;
        }
        osRef.flush();

        scalar totalQ = gSum(Q_.primitiveFieldRef()*mesh_.V());
        Info << "Total energy absorbed [W]: " << totalQ << endl;

        if (mesh_.time().outputTime())
        {
             reflectingCellsVol.write();
        }
    }

    // Clear and initialise the cloud
    // NOTE: Possible to reset original particles, but this requires
    // data transfer for the cloud in differet processors.
    initialise();
}


Foam::tmp<Foam::volScalarField> Foam::radiation::laserDTRM::Rp() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "zero",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar
            (
                "zero",
                dimPower/dimVolume/pow4(dimTemperature),
                0.0
            )
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::radiation::laserDTRM::Ru() const
{
    return Q_.internalField();
}


// ************************************************************************* //
