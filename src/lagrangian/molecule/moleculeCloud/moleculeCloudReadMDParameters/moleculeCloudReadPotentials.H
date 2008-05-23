IOdictionary potentialDict
(
    IOobject
    (
        "potentialDict",
        mesh_.time().system(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// ****************************************************************************
// Pair potentials

Info << nl << "Creating pair potentials:" << nl << endl;

if (!potentialDict.found("pair"))
{
    FatalErrorIn("moleculeCloudReadPotentials.H") << nl
            << "pair potential specification subDict not found"
            << abort(FatalError);
}

const dictionary& pairDict = potentialDict.subDict("pair");

pairPotentials_.setSizeAndNIds(nIds());

label a;

label b;

rCutMax_ = 0.0;

for
(
    a = 0;
    a < idList_.size();
    ++a
)
{
    word idA = idList_[a];

    for
    (
        b = a;
        b < idList_.size();
        ++b
    )
    {
        word idB = idList_[b];

        word pairPotentialName;

        if (a==b)
        {


            if(pairDict.found(idA+"-"+idB))
            {
                pairPotentialName = idA+"-"+idB;
            }

            else
            {
                FatalErrorIn("moleculeCloudReadPotentials.H") << nl
                        << "Pair pairPotential specification subDict "
                        << idA+"-"+idB << " not found"
                        << abort(FatalError);
            }
        }

        else
        {
            if(pairDict.found(idA+"-"+idB))
            {
                pairPotentialName = idA+"-"+idB;
            }

            else if(pairDict.found(idB+"-"+idA))
            {
                pairPotentialName = idB+"-"+idA;
            }

            else
            {
                FatalErrorIn("moleculeCloudReadPotentials.H") << nl
                        << "Pair pairPotential specification subDict "
                        << idA+"-"+idB << " or " << idB+"-"+idA << " not found"
                        << abort(FatalError);
            }

            if
            (
                pairDict.found(idA+"-"+idB)
                && pairDict.found(idB+"-"+idA)
            )
            {
                FatalErrorIn("moleculeCloudReadPotentials.H") << nl
                        << "Pair pairPotential specification subDict "
                        << idA+"-"+idB << " and "
                        << idB+"-"+idA << " found - multiple definition"
                        << abort(FatalError);
            }
        }

        word pairPotentialType
        (
            pairDict.subDict(pairPotentialName).lookup("potentialType")
        );

        if
        (
            pairPotentialType == "maitlandSmithTabulated"
        )
        {
            scalar rCut = readScalar
            (
                pairDict.subDict(pairPotentialName).lookup("rCut")
            );

            if(rCut > rCutMax_)
            {
                rCutMax_ = rCut;
            }

            pairPotentials_.addPotential
            (
                a,
                b,
                pairPotential
                (
                    pairDict.subDict(pairPotentialName)
                )
            );
        }
        else
        {
            FatalErrorIn("moleculeCloudReadPotentials.H") << nl
                    << "pairPotentialType "<< pairPotentialType << " unknown"
                    << abort(FatalError);
        }

        Info << pairPotentialName << ": "
            << pairPotentials_.pairPotentialFunction(a,b) << endl;
    }
}

rCutMax_ += guardRadius_;

rCutMaxSqr_ = rCutMax_ * rCutMax_;

Info << nl << "rCutMax + guardRadius = " << rCutMax_ << endl;

if (potentialDict.found("removalOrder"))
{
    List<word> remOrd = potentialDict.lookup("removalOrder");

    removalOrder_.setSize(remOrd.size());

    forAll(removalOrder_, rO)
    {
        removalOrder_[rO] = findIndex(idList_, remOrd[rO]);
    }
}

// ****************************************************************************
// Tether potentials

iterator mol(this->begin());

DynamicList<label> tetherIds;

for
(
 mol = this->begin();
 mol != this->end();
 ++mol
)
{
    if (mol().tethered())
    {
        if (findIndex(tetherIds, mol().id()) == -1)
        {
            tetherIds.append
            (
                mol().id()
            );
        }
    }
}

if (Pstream::parRun())
{
    List< labelList > allTetherIds(Pstream::nProcs());

    allTetherIds[Pstream::myProcNo()] = tetherIds;

    Pstream::gatherList(allTetherIds);

    if (Pstream::master())
    {
        DynamicList<label> globalTetherIds;

        forAll(allTetherIds, procN)
        {
            const labelList& procNTetherIds = allTetherIds[procN];

            forAll(procNTetherIds, id)
            {
                if (findIndex(globalTetherIds, procNTetherIds[id]) == -1)
                {
                    globalTetherIds.append
                    (
                        procNTetherIds[id]
                    );
                }
            }
        }

        globalTetherIds.shrink();

        tetherIds = globalTetherIds;
    }

    Pstream::scatter(tetherIds);
}

tetherIds.shrink();

if (tetherIds.size())
{
    Info << nl << "Creating tether potentials:" << endl;

    if (!potentialDict.found("tether"))
    {
        FatalErrorIn("moleculeCloudReadPotentials.H") << nl
                << "tether potential specification subDict not found"
                << abort(FatalError);
    }

    const dictionary& tetherDict = potentialDict.subDict("tether");

    tetherPotentials_.setSizeAndTetherIds(tetherIds);

    forAll (tetherIds, tid)
    {
        word tetherPotentialName = idList_[tetherIds[tid]];

        if (!tetherDict.found(tetherPotentialName))
        {
            FatalErrorIn("moleculeCloudReadPotentials.H") << nl
                    << "tether potential specification subDict "
                    << tetherPotentialName << " not found"
                    << abort(FatalError);
        }

        Info << nl << tetherPotentialName << endl;

        word tetherPotentialType =
            tetherDict.subDict(tetherPotentialName).lookup("potentialType");

        if
        (
            tetherPotentialType == "harmonicSpring"
        )
        {

            Info << "\tpotentialType = " << tetherPotentialType << endl;

            scalar springConstant
            (
                readScalar
                (
                    tetherDict.subDict(tetherPotentialName).lookup
                    (
                        "springConstant"
                    )
                )
            );

            Info << "\tspringConstant = " << springConstant << endl;

            tetherPotentials_.addPotential
            (
                tetherIds[tid],
                tetherPotential
                (
                    tetherPotentialName,
                    tetherPotentialType,
                    springConstant
                )
            );
        }

        else
        {
            FatalErrorIn("moleculeCloudReadPotentials.H") << nl
                    << "tetherPotentialType "
                    << tetherPotentialType << " unknown"
                    << abort(FatalError);
        }
    }
}

// ****************************************************************************
// External Forces

gravity_ = vector::zero;

if (potentialDict.found("external"))
{

    Info << nl << "Reading external forces:" << endl;

    const dictionary& externalDict = potentialDict.subDict("external");

    // ************************************************************************
    // gravity

    if (externalDict.found("gravity"))
    {
        gravity_ = externalDict.lookup("gravity");
    }
}

Info << nl << tab << "gravity = " << gravity_ << endl;
