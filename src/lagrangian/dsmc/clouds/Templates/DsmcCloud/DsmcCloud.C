/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "DsmcCloud.H"
#include "BinaryCollisionModel.H"
#include "WallInteractionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::scalar Foam::DsmcCloud<ParcelType>::kb = 1.380650277e-23;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::buildConstProps()
{
    Info<< nl << "typeIds found " << typeIdList_ << endl;

    constProps_.setSize(typeIdList_.size());

    dictionary moleculeProperties
    (
        particleProperties_.subDict("moleculeProperties")
    );

    forAll(typeIdList_, i)
    {
        const word& id(typeIdList_[i]);

        const dictionary& molDict(moleculeProperties.subDict(id));

        constProps_[i] =
        typename ParcelType::constantProperties::constantProperties(molDict);
    }
}


template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::buildCellOccupancy()
{

    forAll(cellOccupancy_, cO)
    {
        cellOccupancy_[cO].clear();
    }

    forAllIter(typename DsmcCloud<ParcelType>, *this, iter)
    {
        cellOccupancy_[iter().cell()].append(&iter());
    }
}


template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::initialise
(
    const IOdictionary& dsmcInitialiseDict
)
{
    const scalar temperature
    (
        readScalar(dsmcInitialiseDict.lookup("temperature"))
    );

    const vector velocity(dsmcInitialiseDict.lookup("velocity"));

    List<word> molecules
    (
        dsmcInitialiseDict.lookup("molecules")
    );

    Field<scalar> numberDensities
    (
        dsmcInitialiseDict.lookup("numberDensities")
    );

    if(molecules.size() != numberDensities.size())
    {
        FatalErrorIn("Foam::Foam::DsmcCloud<ParcelType>::initialise")
            << "molecules and numberDensities must be the same size."
            << nl << abort(FatalError);
    }

    numberDensities /= nParticle_;

    scalar x0 = mesh_.bounds().min().x();
    scalar xR = mesh_.bounds().max().x() - x0;

    scalar y0 = mesh_.bounds().min().y();
    scalar yR = mesh_.bounds().max().y() - y0;

    scalar z0 = mesh_.bounds().min().z();
    scalar zR = mesh_.bounds().max().z() - z0;

    forAll(molecules, i)
    {
        const word& moleculeName(molecules[i]);

        label typeId(findIndex(typeIdList_, moleculeName));

        if(typeId == -1)
        {
            FatalErrorIn("Foam::DsmcCloud<ParcelType>::initialise")
                << "typeId " << moleculeName << "not defined." << nl
                << abort(FatalError);
        }

        const typename ParcelType::constantProperties& cP = constProps(typeId);

        scalar numberDensity = numberDensities[i];

        scalar spacing = pow(numberDensity,-(1.0/3.0));

        int ni = label(xR/spacing) + 1;
        int nj = label(yR/spacing) + 1;
        int nk = label(zR/spacing) + 1;

        vector delta(xR/ni, yR/nj, zR/nk);

        scalar pert = spacing;

        for (int i = 0; i < ni; i++)
        {
            for (int j = 0; j < nj; j++)
            {
                for (int k = 0; k < nk; k++)
                {
                    point p
                    (
                        x0 + (i + 0.5)*delta.x(),
                        y0 + (j + 0.5)*delta.y(),
                        z0 + (k + 0.5)*delta.z()
                    );

                    p.x() += pert*(rndGen_.scalar01() - 0.5);
                    p.y() += pert*(rndGen_.scalar01() - 0.5);
                    p.z() += pert*(rndGen_.scalar01() - 0.5);

                    label cell = mesh_.findCell(p);

                    vector U = equipartitionLinearVelocity
                    (
                        temperature,
                        cP.mass()
                    );

                    U += velocity;

                    if (cell >= 0)
                    {
                        addNewParcel
                        (
                            p,
                            U,
                            cell,
                            typeId
                        );
                    }
                }
            }
        }
    }
}


template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::collisions()
{
    Info<< "DsmcCloud collisions() - PLACEHOLDER" << endl;
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::addNewParcel
(
    const vector& position,
    const vector& U,
    const label cellId,
    const label typeId
)
{
    ParcelType* pPtr = new ParcelType
    (
        *this,
        position,
        U,
        cellId,
        typeId
    );

    addParticle(pPtr);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::DsmcCloud<ParcelType>::DsmcCloud
(
    const word& cloudType,
    const fvMesh& mesh
)
:
    Cloud<ParcelType>(mesh, cloudType, false),
    DsmcBaseCloud(),
    cloudType_(cloudType),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            cloudType + "Properties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    nParticle_(readScalar(particleProperties_.lookup("nEquivalentParticles"))),
    cellOccupancy_(mesh_.nCells()),
    constProps_(),
    rndGen_(label(971501)),
    binaryCollisionModel_
    (
        BinaryCollisionModel<DsmcCloud<ParcelType> >::New
        (
            particleProperties_,
            *this
        )
    ),
    wallInteractionModel_
    (
        WallInteractionModel<DsmcCloud<ParcelType> >::New
        (
            particleProperties_,
            *this
        )
    )
{
    buildConstProps();

    buildCellOccupancy();
}


template<class ParcelType>
Foam::DsmcCloud<ParcelType>::DsmcCloud
(
    const word& cloudType,
    const fvMesh& mesh,
    const IOdictionary& dsmcInitialiseDict
)
    :
    Cloud<ParcelType>(mesh, cloudType, false),
    DsmcBaseCloud(),
    cloudType_(cloudType),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            cloudType + "Properties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    nParticle_(readScalar(particleProperties_.lookup("nEquivalentParticles"))),
    cellOccupancy_(),
    constProps_(),
    rndGen_(label(971501)),
    binaryCollisionModel_(),
    wallInteractionModel_()
{
    clear();

    buildConstProps();

    initialise(dsmcInitialiseDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::DsmcCloud<ParcelType>::~DsmcCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::evolve()
{
    typename ParcelType::trackData td(*this);

    //this->injection().inject(td);

    if (debug)
    {
        this->dumpParticlePositions();
    }

    // Move the particles ballistically with their current velocities
    Cloud<ParcelType>::move(td);

    // Calculate new velocities via stochastic collisions
    collisions();
}


template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::info() const
{
    Info<< "Cloud name: " << this->name() << nl
        << "    Current number of parcels       = "
        << returnReduce(this->size(), sumOp<label>()) << nl
        << "    Current mass in system          = "
        << returnReduce(massInSystem(), sumOp<scalar>()) << nl
        << endl;
}


template<class ParcelType>
Foam::vector Foam::DsmcCloud<ParcelType>::equipartitionLinearVelocity
(
    scalar temperature,
    scalar mass
)
{
    return sqrt(kb*temperature/mass)*vector
    (
        rndGen_.GaussNormal(),
        rndGen_.GaussNormal(),
        rndGen_.GaussNormal()
    );
}


template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::dumpParticlePositions() const
{
    OFstream pObj
    (
        this->db().time().path()/"parcelPositions_"
      + this->name() + "_"
      + this->db().time().timeName() + ".obj"
    );

    forAllConstIter(typename DsmcCloud<ParcelType>, *this, iter)
    {
        const ParcelType& p = iter();

        pObj<< "v " << p.position().x()
            << " "  << p.position().y()
            << " "  << p.position().z()
            << nl;
    }

    pObj.flush();
}


// ************************************************************************* //
