/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "KinematicSurfaceFilm.H"
#include "surfaceFilmRegionModel.H"
#include "liquidFilmModel.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "Pstream.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CloudType>
Foam::wordList Foam::KinematicSurfaceFilm<CloudType>::interactionTypeNames_
{
    "absorb", "bounce", "splashBai"
};


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class CloudType>
typename Foam::KinematicSurfaceFilm<CloudType>::interactionType
Foam::KinematicSurfaceFilm<CloudType>::interactionTypeEnum(const word& it) const
{
    forAll(interactionTypeNames_, i)
    {
        if (interactionTypeNames_[i] == it)
        {
            return interactionType(i);
        }
    }

    FatalErrorInFunction
        << "Unknown interaction type " << it
        << ". Valid interaction types include: " << interactionTypeNames_
        << abort(FatalError);

    return interactionType(0);
}


template<class CloudType>
Foam::word Foam::KinematicSurfaceFilm<CloudType>::interactionTypeStr
(
    const interactionType& it
) const
{
    if (it >= interactionTypeNames_.size())
    {
        FatalErrorInFunction
            << "Unknown interaction type enumeration" << abort(FatalError);
    }

    return interactionTypeNames_[it];
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::vector Foam::KinematicSurfaceFilm<CloudType>::tangentVector
(
    const vector& v
) const
{
    vector tangent(Zero);
    scalar magTangent = 0.0;

    while (magTangent < SMALL)
    {
        const vector vTest(rndGen_.sample01<vector>());
        tangent = vTest - (vTest & v)*v;
        magTangent = mag(tangent);
    }

    return tangent/magTangent;
}


template<class CloudType>
Foam::vector Foam::KinematicSurfaceFilm<CloudType>::splashDirection
(
    const vector& tanVec1,
    const vector& tanVec2,
    const vector& nf
) const
{
    // Azimuthal angle [rad]
    const scalar phiSi = twoPi*rndGen_.sample01<scalar>();

    // Ejection angle [rad]
    const scalar thetaSi = degToRad(rndGen_.sample01<scalar>()*(50 - 5) + 5);

    // Direction vector of new parcel
    const scalar alpha = sin(thetaSi);
    const scalar dcorr = cos(thetaSi);
    const vector normal(alpha*(tanVec1*cos(phiSi) + tanVec2*sin(phiSi)));
    vector dirVec(dcorr*nf);
    dirVec += normal;

    return dirVec/mag(dirVec);
}


template<class CloudType>
void Foam::KinematicSurfaceFilm<CloudType>::initFilmModels()
{
    const fvMesh& mesh = this->owner().mesh();

    // set up filmModel pointer
    if (!filmModel_)
    {
        filmModel_ =
            const_cast<regionFilm*>
            (
                mesh.time().objectRegistry::template findObject
                <
                    regionFilm
                >
                (
                    "surfaceFilmProperties"
                )
            );
    }

    if (areaFilms_.size() == 0)
    {
        // set up areaFilms
        const wordList names =
            mesh.time().objectRegistry::template
                sortedNames<regionModels::regionFaModel>();

        forAll(names, i)
        {
            const regionModels::regionFaModel* regionFa =
                mesh.time().objectRegistry::template findObject
                <
                    regionModels::regionFaModel
                >(names[i]);

            if (regionFa && isA<areaFilm>(*regionFa))
            {
                areaFilm& film =
                    const_cast<areaFilm&>(refCast<const areaFilm>(*regionFa));
                areaFilms_.append(&film);
            }
        }
    }
}


template<class CloudType>
void Foam::KinematicSurfaceFilm<CloudType>::init(bool binitThermo)
{
    if (binitThermo)
    {
        this->coeffDict().readEntry("pRef", pRef_);
        this->coeffDict().readEntry("TRef", TRef_);
        thermo_ = new liquidMixtureProperties(this->coeffDict().subDict("thermo"));
    }
}


template<class CloudType>
template<class filmType>
void Foam::KinematicSurfaceFilm<CloudType>::absorbInteraction
(
    filmType& film,
    const parcelType& p,
    const polyPatch& pp,
    const label facei,
    const scalar mass,
    bool& keepParticle
)
{
    if (debug)
    {
        Info<< "Parcel " << p.origId() << " absorbInteraction" << endl;
    }

    // Patch face normal
    const vector& nf = pp.faceNormals()[facei];

    // Patch velocity
    const vector& Up = this->owner().U().boundaryField()[pp.index()][facei];

    // Relative parcel velocity
    const vector Urel(p.U() - Up);

    // Parcel normal velocity
    const vector Un(nf*(Urel & nf));

    // Parcel tangential velocity
    const vector Ut(Urel - Un);

    film.addSources
    (
        pp.index(),
        facei,
        mass,                             // mass
        mass*Ut,                          // tangential momentum
        mass*mag(Un),                     // impingement pressure
        0                                 // energy
    );

    this->nParcelsTransferred()++;

    this->totalMassTransferred() += mass;

    keepParticle = false;
}


template<class CloudType>
void Foam::KinematicSurfaceFilm<CloudType>::bounceInteraction
(
    parcelType& p,
    const polyPatch& pp,
    const label facei,
    bool& keepParticle
) const
{
    if (debug)
    {
        Info<< "Parcel " << p.origId() << " bounceInteraction" << endl;
    }

    // Patch face normal
    const vector& nf = pp.faceNormals()[facei];

    // Patch velocity
    const vector& Up = this->owner().U().boundaryField()[pp.index()][facei];

    // Relative parcel velocity
    const vector Urel(p.U() - Up);

    // Flip parcel normal velocity component
    p.U() -= 2.0*nf*(Urel & nf);

    keepParticle = true;
}


template<class CloudType>
template<class filmType>
void Foam::KinematicSurfaceFilm<CloudType>::drySplashInteraction
(
    filmType& filmModel,
    const scalar sigma,
    const scalar mu,
    const parcelType& p,
    const polyPatch& pp,
    const label facei,
    bool& keepParticle
)
{
    if (debug)
    {
        Info<< "Parcel " << p.origId() << " drySplashInteraction" << endl;
    }

    // Patch face velocity and normal
    const vector& Up = this->owner().U().boundaryField()[pp.index()][facei];
    const vector& nf = pp.faceNormals()[facei];

    // Local pressure
    //const scalar pc = thermo_.thermo().p()[p.cell()];

    // Retrieve parcel properties
    const scalar m = p.mass()*p.nParticle();
    const scalar rho = p.rho();
    const scalar d = p.d();
    const vector Urel(p.U() - Up);
    const vector Un(nf*(Urel & nf));

    // Laplace number
    const scalar La = rho*sigma*d/sqr(mu);

    // Weber number
    const scalar We = rho*magSqr(Un)*d/sigma;

    // Critical Weber number
    const scalar Wec = Adry_*pow(La, -0.183);

    if (We < Wec) // Adhesion - assume absorb
    {
        absorbInteraction<filmType>
            (filmModel, p, pp, facei, m, keepParticle);
    }
    else // Splash
    {
        // Ratio of incident mass to splashing mass
        const scalar mRatio = 0.2 + 0.6*rndGen_.sample01<scalar>();
        splashInteraction<filmType>
            (filmModel, p, pp, facei, mRatio, We, Wec, sigma, keepParticle);
    }
}


template<class CloudType>
template<class filmType>
void Foam::KinematicSurfaceFilm<CloudType>::wetSplashInteraction
(
    filmType& filmModel,
    const scalar sigma,
    const scalar mu,
    parcelType& p,
    const polyPatch& pp,
    const label facei,
    bool& keepParticle
)
{
    if (debug)
    {
        Info<< "Parcel " << p.origId() << " wetSplashInteraction" << endl;
    }

    // Patch face velocity and normal
    const vector& Up = this->owner().U().boundaryField()[pp.index()][facei];
    const vector& nf = pp.faceNormals()[facei];

    // Retrieve parcel properties
    const scalar m = p.mass()*p.nParticle();
    const scalar rho = p.rho();
    const scalar d = p.d();
    vector& U = p.U();
    const vector Urel(p.U() - Up);
    const vector Un(nf*(Urel & nf));
    const vector Ut(Urel - Un);

    // Laplace number
    const scalar La = rho*sigma*d/sqr(mu);

    // Weber number
    const scalar We = rho*magSqr(Un)*d/sigma;

    // Critical Weber number
    const scalar Wec = Awet_*pow(La, -0.183);

    if (We < 2) // Adhesion - assume absorb
    {
        absorbInteraction<filmType>
            (filmModel, p, pp, facei, m, keepParticle);
    }
    else if ((We >= 2) && (We < 20)) // Bounce
    {
        // Incident angle of impingement
        const scalar theta = piByTwo - acos(U/mag(U) & nf);

        // Restitution coefficient
        const scalar epsilon = 0.993 - theta*(1.76 - theta*(1.56 - theta*0.49));

        // Update parcel velocity
        U = -epsilon*(Un) + 5.0/7.0*(Ut);

        keepParticle = true;
        return;
    }
    else if ((We >= 20) && (We < Wec)) // Spread - assume absorb
    {
        absorbInteraction<filmType>
            (filmModel, p, pp, facei, m, keepParticle);
    }
    else    // Splash
    {
        // Ratio of incident mass to splashing mass
        // splash mass can be > incident mass due to film entrainment
        const scalar mRatio = 0.2 + 0.9*rndGen_.sample01<scalar>();
        splashInteraction<filmType>
            (filmModel, p, pp, facei, mRatio, We, Wec, sigma, keepParticle);
    }
}


template<class CloudType>
template<class filmType>
void Foam::KinematicSurfaceFilm<CloudType>::splashInteraction
(
    filmType& filmModel,
    const parcelType& p,
    const polyPatch& pp,
    const label facei,
    const scalar mRatio,
    const scalar We,
    const scalar Wec,
    const scalar sigma,
    bool& keepParticle
)
{
    // Patch face velocity and normal
    const fvMesh& mesh = this->owner().mesh();
    const vector& Up = this->owner().U().boundaryField()[pp.index()][facei];
    const vector& nf = pp.faceNormals()[facei];

    // Determine direction vectors tangential to patch normal
    const vector tanVec1(tangentVector(nf));
    const vector tanVec2(nf^tanVec1);

    // Retrieve parcel properties
    const scalar np = p.nParticle();
    const scalar m = p.mass()*np;
    const scalar d = p.d();
    const vector Urel(p.U() - Up);
    const vector Un(nf*(Urel & nf));
    const vector Ut(Urel - Un);
    const vector& posC = mesh.C()[p.cell()];
    const vector& posCf = mesh.Cf().boundaryField()[pp.index()][facei];

    // Total mass of (all) splashed parcels
    const scalar mSplash = m*mRatio;

    // Number of splashed particles per incoming particle
    const scalar Ns = 5.0*(We/Wec - 1.0);

    // Average diameter of splashed particles
    const scalar dBarSplash = 1/cbrt(6.0)*cbrt(mRatio/Ns)*d + ROOTVSMALL;

    // Cumulative diameter splash distribution
    const scalar dMax = 0.9*cbrt(mRatio)*d;
    const scalar dMin = 0.1*dMax;
    const scalar K = exp(-dMin/dBarSplash) - exp(-dMax/dBarSplash);

    // Surface energy of secondary parcels [J]
    scalar ESigmaSec = 0;

    // Sample splash distribution to determine secondary parcel diameters
    scalarList dNew(parcelsPerSplash_);
    scalarList npNew(parcelsPerSplash_);
    forAll(dNew, i)
    {
        const scalar y = rndGen_.sample01<scalar>();
        dNew[i] = -dBarSplash*log(exp(-dMin/dBarSplash) - y*K);
        npNew[i] = mRatio*np*pow3(d)/pow3(dNew[i])/parcelsPerSplash_;
        ESigmaSec += npNew[i]*sigma*p.areaS(dNew[i]);
    }

    // Incident kinetic energy [J]
    const scalar EKIn = 0.5*m*magSqr(Un);

    // Incident surface energy [J]
    const scalar ESigmaIn = np*sigma*p.areaS(d);

    // Dissipative energy
    const scalar Ed = max(0.8*EKIn, np*Wec/12*pi*sigma*sqr(d));

    // Total energy [J]
    const scalar EKs = EKIn + ESigmaIn - ESigmaSec - Ed;

    // Switch to absorb if insufficient energy for splash
    if (EKs <= 0)
    {
        absorbInteraction<filmType>
            (filmModel, p, pp, facei, m, keepParticle);
        return;
    }

    // Helper variables to calculate magUns0
    const scalar logD = log(d);
    const scalar coeff2 = log(dNew[0]) - logD + ROOTVSMALL;
    scalar coeff1 = 0.0;
    forAll(dNew, i)
    {
        coeff1 += sqr(log(dNew[i]) - logD);
    }

    // Magnitude of the normal velocity of the first splashed parcel
    const scalar magUns0 =
        sqrt(2.0*parcelsPerSplash_*EKs/mSplash/(1.0 + coeff1/sqr(coeff2)));

    // Set splashed parcel properties
    forAll(dNew, i)
    {
        const vector dirVec = splashDirection(tanVec1, tanVec2, -nf);

        // Create a new parcel by copying source parcel
        parcelType* pPtr = new parcelType(p);

        pPtr->origId() = pPtr->getNewParticleID();

        pPtr->origProc() = Pstream::myProcNo();

        if (splashParcelType_ >= 0)
        {
            pPtr->typeId() = splashParcelType_;
        }

        // Perturb new parcels towards the owner cell centre
        pPtr->track(0.5*rndGen_.sample01<scalar>()*(posC - posCf), 0);

        pPtr->nParticle() = npNew[i];

        pPtr->d() = dNew[i];

        pPtr->U() = dirVec*(mag(Cf_*Ut) + magUns0*(log(dNew[i]) - logD)/coeff2);

        // Apply correction to velocity for 2-D cases
        meshTools::constrainDirection(mesh, mesh.solutionD(), pPtr->U());

        // Add the new parcel
        this->owner().addParticle(pPtr);

        nParcelsSplashed_++;
    }

    // Transfer remaining part of parcel to film 0 - splashMass can be -ve
    // if entraining from the film
    const scalar mDash = m - mSplash;
    absorbInteraction<filmType>
        (filmModel, p, pp, facei, mDash, keepParticle);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::KinematicSurfaceFilm<CloudType>::KinematicSurfaceFilm
(
    const dictionary& dict,
    CloudType& owner,
    const word& type,
    bool initThermo
)
:
    SurfaceFilmModel<CloudType>(dict, owner, type),
    rndGen_(owner.rndGen()),
    thermo_(nullptr),
    filmModel_(nullptr),
    areaFilms_(0),
    interactionType_
    (
        interactionTypeEnum(this->coeffDict().getWord("interactionType"))
    ),
    deltaWet_(0.0),
    splashParcelType_(0),
    parcelsPerSplash_(0),
    Adry_(0.0),
    Awet_(0.0),
    Cf_(0.0),
    nParcelsSplashed_(0)
{
    Info<< "    Applying " << interactionTypeStr(interactionType_)
        << " interaction model" << endl;

    if (interactionType_ == itSplashBai)
    {
        this->coeffDict().readEntry("deltaWet", deltaWet_);
        splashParcelType_ =
            this->coeffDict().getOrDefault("splashParcelType", -1);
        parcelsPerSplash_ =
            this->coeffDict().getOrDefault("parcelsPerSplash", 2);
        this->coeffDict().readEntry("Adry", Adry_);
        this->coeffDict().readEntry("Awet", Awet_);
        this->coeffDict().readEntry("Cf", Cf_);
        init(initThermo);
    }
}


template<class CloudType>
Foam::KinematicSurfaceFilm<CloudType>::KinematicSurfaceFilm
(
    const KinematicSurfaceFilm<CloudType>& sfm,
    bool initThermo
)
:
    SurfaceFilmModel<CloudType>(sfm),
    rndGen_(sfm.rndGen_),
    thermo_(nullptr),
    filmModel_(nullptr),
    areaFilms_(0),
    interactionType_(sfm.interactionType_),
    deltaWet_(sfm.deltaWet_),
    splashParcelType_(sfm.splashParcelType_),
    parcelsPerSplash_(sfm.parcelsPerSplash_),
    Adry_(sfm.Adry_),
    Awet_(sfm.Awet_),
    Cf_(sfm.Cf_),
    nParcelsSplashed_(sfm.nParcelsSplashed_)
{
    if (interactionType_ == itSplashBai)
    {
        init(initThermo);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::KinematicSurfaceFilm<CloudType>::transferParcel
(
    parcelType& p,
    const polyPatch& pp,
    bool& keepParticle
)
{
    const label patchi = pp.index();

    bool bInteraction(false);

    initFilmModels();

    // Check the singleLayer film models
    if (filmModel_)
    {
        if (filmModel_->isRegionPatch(patchi))
        {
            const label facei = pp.whichFace(p.face());

            switch (interactionType_)
            {
                case itBounce:
                {
                    bounceInteraction(p, pp, facei, keepParticle);

                    break;
                }
                case itAbsorb:
                {
                    const scalar m = p.nParticle()*p.mass();

                    absorbInteraction<regionFilm>
                        (*filmModel_, p, pp, facei, m, keepParticle);

                    break;
                }
                case itSplashBai:
                {
                    bool dry = this->deltaFilmPatch_[patchi][facei] < deltaWet_;

                    const scalarField X(thermo_->size(), 1);
                    const scalar mu = thermo_->mu(pRef_, TRef_, X);
                    const scalar sigma = thermo_->sigma(pRef_, TRef_, X);

                    if (dry)
                    {
                        drySplashInteraction<regionFilm>
                            (*filmModel_, sigma, mu, p, pp, facei, keepParticle);
                    }
                    else
                    {
                        wetSplashInteraction<regionFilm>
                            (*filmModel_, sigma, mu, p, pp, facei, keepParticle);
                    }

                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Unknown interaction type enumeration"
                        << abort(FatalError);
                }
            }

            // Transfer parcel/parcel interactions complete
            bInteraction = true;
        }
    }


    for (areaFilm& film : areaFilms_)
    {
        if (patchi == film.patchID())
        {
            const label facei = pp.whichFace(p.face());

            switch (interactionType_)
            {
                // It only supports absorp model
                case itAbsorb:
                {
                    const scalar m = p.nParticle()*p.mass();

                    absorbInteraction<areaFilm>
                    (
                        film, p, pp, facei, m, keepParticle
                    );
                    break;
                }
                case itBounce:
                {
                    bounceInteraction(p, pp, facei, keepParticle);

                    break;
                }
                case itSplashBai:
                {
                    bool dry = film.h()[facei] < deltaWet_;

                    regionModels::areaSurfaceFilmModels::liquidFilmModel& liqFilm =
                        refCast
                        <   regionModels::areaSurfaceFilmModels::liquidFilmModel
                        >(film);

                    const scalarField X(liqFilm.thermo().size(), 1);
                    const scalar pRef = film.pRef();
                    const scalar TRef = liqFilm.Tref();

                    const scalar mu = liqFilm.thermo().mu(pRef, TRef, X);
                    const scalar sigma =
                        liqFilm.thermo().sigma(pRef, TRef, X);

                    if (dry)
                    {
                        drySplashInteraction<areaFilm>
                            (film, sigma, mu, p, pp, facei, keepParticle);
                    }
                    else
                    {
                        wetSplashInteraction<areaFilm>
                            (film, sigma, mu, p, pp, facei, keepParticle);
                    }

                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Unknown interaction type enumeration"
                        << abort(FatalError);
                }
            }
            // Transfer parcel/parcel interactions complete
            bInteraction = true;
        }
    }

    // Parcel not interacting with film
    return bInteraction;
}


template<class CloudType>
void Foam::KinematicSurfaceFilm<CloudType>::cacheFilmFields
(
    const label filmPatchi,
    const label primaryPatchi,
    const regionModels::surfaceFilmModels::surfaceFilmRegionModel& filmModel
)
{
    SurfaceFilmModel<CloudType>::cacheFilmFields
    (
        filmPatchi,
        primaryPatchi,
        filmModel
    );
}


template<class CloudType>
void Foam::KinematicSurfaceFilm<CloudType>::cacheFilmFields
(
    const label filmPatchi,
    const areaFilm& filmModel
)
{
    SurfaceFilmModel<CloudType>::cacheFilmFields
    (
        filmPatchi,
        filmModel
    );
}


template<class CloudType>
void Foam::KinematicSurfaceFilm<CloudType>::setParcelProperties
(
    parcelType& p,
    const label filmFacei
) const
{
    SurfaceFilmModel<CloudType>::setParcelProperties(p, filmFacei);
}


template<class CloudType>
void Foam::KinematicSurfaceFilm<CloudType>::info(Ostream& os)
{
    SurfaceFilmModel<CloudType>::info(os);

    label nSplash0 = this->template getModelProperty<label>("nParcelsSplashed");
    label nSplashTotal =
        nSplash0 + returnReduce(nParcelsSplashed_, sumOp<label>());

    os  << "      - new splash parcels          = " << nSplashTotal << endl;

    if (this->writeTime())
    {
        this->setModelProperty("nParcelsSplashed", nSplashTotal);
        nParcelsSplashed_ = 0;
    }
}


// ************************************************************************* //
