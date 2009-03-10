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
#include "InflowBoundaryModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::scalar Foam::DsmcCloud<ParcelType>::kb = 1.380650277e-23;

template<class ParcelType>
Foam::scalar Foam::DsmcCloud<ParcelType>::Tref = 273;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::buildConstProps()
{
    Info<< nl << "Constructing constant properties for" << endl;
    constProps_.setSize(typeIdList_.size());

    dictionary moleculeProperties
    (
        particleProperties_.subDict("moleculeProperties")
    );

    forAll(typeIdList_, i)
    {
        const word& id(typeIdList_[i]);

        Info<< "    " << id << endl;

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
    Info<< nl << "Initialising particles" << endl;

    const scalar temperature
    (
        readScalar(dsmcInitialiseDict.lookup("temperature"))
    );

    const vector velocity(dsmcInitialiseDict.lookup("velocity"));

    const dictionary& numberDensitiesDict
    (
        dsmcInitialiseDict.subDict("numberDensities")
    );

    List<word> molecules(numberDensitiesDict.toc());

    Field<scalar> numberDensities(molecules.size());

    forAll(molecules, i)
    {
        numberDensities[i] = readScalar
        (
            numberDensitiesDict.lookup(molecules[i])
        );
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

                    scalar Ei = equipartitionInternalEnergy
                    (
                        temperature,
                        cP.internalDegreesOfFreedom()
                    );

                    U += velocity;

                    if (cell >= 0)
                    {
                        addNewParcel
                        (
                            p,
                            U,
                            Ei,
                            cell,
                            typeId
                        );
                    }
                }
            }
        }
    }

    // Initialise the sigmaTcRMax_ field to the product of the cross section of
    // the most abundant species and the most probable thermal speed (Bird,
    // p222-223)

    label mostAbundantType(findMax(numberDensities));

    const typename ParcelType::constantProperties& cP = constProps
    (
        mostAbundantType
    );

    sigmaTcRMax_.internalField() = cP.sigmaT()*maxwellianMostProbableSpeed
    (
        temperature,
        cP.mass()
    );

    sigmaTcRMax_.correctBoundaryConditions();
}


template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::collisions()
{
    buildCellOccupancy();

    // Temporary storage for subCells
    List<DynamicList<label> > subCells(8);

    scalar deltaT = mesh_.time().deltaT().value();

    label collisionCandidates = 0;

    label collisions = 0;

    forAll(cellOccupancy_, celli)
    {
        const DynamicList<ParcelType*>& cellParcels(cellOccupancy_[celli]);

        label nC(cellParcels.size());

        if (nC > 1)
        {

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // Assign particles to one of 8 Cartesian subCells

            // Clear temporary lists
            forAll(subCells, i)
            {
                subCells[i].clear();
            }

            // Inverse addressing specifying which subCell a parcel is in
            List<label> whichSubCell(cellParcels.size());

            const point& cC = mesh_.cellCentres()[celli];

            forAll(cellParcels, i)
            {
                ParcelType* p = cellParcels[i];

                vector relPos = p->position() - cC;

                label subCell =
                    pos(relPos.x()) + 2*pos(relPos.y()) + 4*pos(relPos.z());

                subCells[subCell].append(i);

                whichSubCell[i] = subCell;
            }

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            scalar sigmaTcRMax = sigmaTcRMax_[celli];

            scalar selectedPairs = collisionSelectionRemainder_[celli]
              + 0.5*nC*(nC-1)*nParticle_*sigmaTcRMax*deltaT
               /mesh_.cellVolumes()[celli];

            label nCandidates(selectedPairs);

            collisionSelectionRemainder_[celli] = selectedPairs - nCandidates;

            collisionCandidates += nCandidates;

            for(label c = 0; c < nCandidates; c++)
            {
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // subCell candidate selection procedure

                // Select the first collision candidate
                label candidateP = rndGen_.integer(0, nC-1);

                // Declare the second collision candidate
                label candidateQ = -1;

                List<label> subCellPs = subCells[whichSubCell[candidateP]];

                label nSC = subCellPs.size();

                if (nSC > 1)
                {
                    // If there are two or more particle in a subCell, choose
                    // another from the same cell.  If the same candidate is
                    // chosen, choose again.

                    do
                    {
                        candidateQ = subCellPs[rndGen_.integer(0, nSC-1)];

                    } while(candidateP == candidateQ);
                }
                else
                {
                    // Select a possible second collision candidate from the
                    // whole cell.  If the same candidate is chosen, choose
                    // again.

                    do
                    {
                        candidateQ = rndGen_.integer(0, nC-1);

                    } while(candidateP == candidateQ);
                }

                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // uniform candidate selection procedure

                // // Select the first collision candidate
                // label candidateP = rndGen_.integer(0, nC-1);

                // // Select a possible second collision candidate
                // label candidateQ = rndGen_.integer(0, nC-1);

                // // If the same candidate is chosen, choose again
                // while(candidateP == candidateQ)
                // {
                //     candidateQ = rndGen_.integer(0, nC-1);
                // }

                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                ParcelType* parcelP = cellParcels[candidateP];
                ParcelType* parcelQ = cellParcels[candidateQ];

                label typeIdP = parcelP->typeId();
                label typeIdQ = parcelQ->typeId();

                scalar sigmaTcR = binaryCollision().sigmaTcR
                (
                    typeIdP,
                    typeIdQ,
                    parcelP->U(),
                    parcelQ->U()
                );

                // Update the maximum value of sigmaTcR stored, but use the
                // initial value in the acceptance-rejection criteria because
                // the number of collision candidates selected was based on this

                if (sigmaTcR > sigmaTcRMax_[celli])
                {
                    sigmaTcRMax_[celli] = sigmaTcR;
                }

                if ((sigmaTcR/sigmaTcRMax) > rndGen_.scalar01())
                {
                    binaryCollision().collide
                    (
                        typeIdP,
                        typeIdQ,
                        parcelP->U(),
                        parcelQ->U(),
                        parcelP->Ei(),
                        parcelQ->Ei()
                    );

                    collisions++;
                }
            }
        }
    }

    reduce(collisions, sumOp<label>());

    reduce(collisionCandidates, sumOp<label>());

    if (collisionCandidates)
    {
        Info<< "    Collisions                      = "
            << collisions << nl
            << "    Acceptance rate                 = "
            << scalar(collisions)/scalar(collisionCandidates) << nl
            << endl;
    }
    else
    {
        Info<< "    No collisions" << endl;
    }
}


template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::resetSurfaceDataFields()
{
    volScalarField::GeometricBoundaryField& qBF(q_.boundaryField());

    forAll(qBF, p)
    {
        qBF[p] = 0.0;
    }

    volVectorField::GeometricBoundaryField& fDBF(fD_.boundaryField());

    forAll(fDBF, p)
    {
        fDBF[p] = vector::zero;
    }
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::addNewParcel
(
    const vector& position,
    const vector& U,
    const scalar Ei,
    const label cellId,
    const label typeId
)
{
    ParcelType* pPtr = new ParcelType
    (
        *this,
        position,
        U,
        Ei,
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
    const volScalarField& T,
    const volVectorField& U
)
:
    Cloud<ParcelType>(T.mesh(), cloudType, false),
    DsmcBaseCloud(),
    cloudType_(cloudType),
    mesh_(T.mesh()),
    particleProperties_
    (
        IOobject
        (
            cloudType + "Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    nParticle_(readScalar(particleProperties_.lookup("nEquivalentParticles"))),
    cellOccupancy_(mesh_.nCells()),
    sigmaTcRMax_
    (
        IOobject
        (
            this->name() + "SigmaTcRMax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    collisionSelectionRemainder_(mesh_.nCells(), 0),
    q_
    (
        IOobject
        (
            this->name() + "q_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, 0, -3, 0, 0), 0.0)
    ),
    fD_
    (
        IOobject
        (
            this->name() + "fD_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, -1, -2, 0, 0),
            vector::zero
        )
    ),
    constProps_(),
    rndGen_(label(149382906) + 7183*Pstream::myProcNo()),
    T_(T),
    U_(U),
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
    ),
    inflowBoundaryModel_
    (
        InflowBoundaryModel<DsmcCloud<ParcelType> >::New
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
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    nParticle_(readScalar(particleProperties_.lookup("nEquivalentParticles"))),
    cellOccupancy_(),
    sigmaTcRMax_
    (
        IOobject
        (
            this->name() + "SigmaTcRMax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, 3, -1, 0, 0), 0.0)
    ),
    collisionSelectionRemainder_(),
    q_
    (
        IOobject
        (
            this->name() + "q_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, 0, -3, 0, 0), 0.0)
    ),
    fD_
    (
        IOobject
        (
            this->name() + "fD_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, -1, -2, 0, 0),
            vector::zero
        )
    ),
    constProps_(),
    rndGen_(label(971501) + 1526*Pstream::myProcNo()),
    T_
    (
        volScalarField
        (
            IOobject
            (
                "T",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero",  dimensionSet(0, 0, 0, 1, 0), 0.0)
        )
    ),
    U_
    (
        volVectorField
        (
            IOobject
            (
                "U",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector
            (
                "zero",
                dimensionSet(0, 1, -1, 0, 0),
                vector::zero
            )
        )
    ),
    binaryCollisionModel_(),
    wallInteractionModel_(),
    inflowBoundaryModel_()
{
    clear();

    buildConstProps();

    IOdictionary dsmcInitialiseDict
    (
        IOobject
        (
            "dsmcInitialiseDict",
            mesh_.time().system(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

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

    // Reset the surface data collection fields
    resetSurfaceDataFields();

    if (debug)
    {
       this->dumpParticlePositions();
    }

    // Insert new particles from the inflow boundary
    this->inflowBoundary().inflow();

    // Move the particles ballistically with their current velocities
    Cloud<ParcelType>::move(td);

    // Calculate new velocities via stochastic collisions
    collisions();
}


template<class ParcelType>
void Foam::DsmcCloud<ParcelType>::info() const
{
    label nDsmcParticles = this->size();
    reduce(nDsmcParticles, sumOp<label>());

    scalar nMol = nDsmcParticles*nParticle_;

    vector linearMomentum = linearMomentumOfSystem();
    reduce(linearMomentum, sumOp<vector>());

    scalar linearKineticEnergy = linearKineticEnergyOfSystem();
    reduce(linearKineticEnergy, sumOp<scalar>());

    scalar internalEnergy = internalEnergyOfSystem();
    reduce(internalEnergy, sumOp<scalar>());

    Info<< "Cloud name: " << this->name() << nl
        << "    Number of dsmc particles        = "
        << nDsmcParticles << nl
        << "    Number of molecules             = "
        << nMol << nl
        << "    Mass in system                  = "
        << returnReduce(massInSystem(), sumOp<scalar>()) << nl
        << "    Average linear momentum         = "
        << linearMomentum/nMol << nl
        << "   |Average linear momentum|        = "
        << mag(linearMomentum)/nMol << nl
        << "    Average linear kinetic energy   = "
        << linearKineticEnergy/nMol << nl
        << "    Average internal energy         = "
        << internalEnergy/nMol << nl
        << "    Average total energy            = "
        << (internalEnergy + linearKineticEnergy)/nMol << nl
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
Foam::scalar Foam::DsmcCloud<ParcelType>::equipartitionInternalEnergy
(
    scalar temperature,
    scalar iDof
)
{
    scalar Ei = 0.0;

    if
    (
        iDof < 2.0 + SMALL
     && iDof > 2.0 - SMALL
    )
    {
        // Special case for iDof = 2, i.e. diatomics;
        Ei = -log(rndGen_.scalar01())*kb*temperature;
    }
    else
    {
        scalar a = 0.5*iDof - 1;

        scalar energyRatio;

        scalar P;

        do
        {
            energyRatio = 10*rndGen_.scalar01();

            P = pow((energyRatio/a), a)*exp(a - energyRatio);

        } while (P < rndGen_.scalar01());

        Ei = energyRatio*kb*temperature;
    }

    return Ei;
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
