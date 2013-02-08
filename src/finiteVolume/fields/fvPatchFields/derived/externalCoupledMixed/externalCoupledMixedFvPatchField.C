/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "externalCoupledMixedFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "IFstream.H"
#include "OFstream.H"
#include "globalIndex.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
Foam::word Foam::externalCoupledMixedFvPatchField<Type>::lockName = "OpenFOAM";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::externalCoupledMixedFvPatchField<Type>::baseDir() const
{
    word regionName(this->dimensionedInternalField().mesh().name());
    if (regionName == polyMesh::defaultRegion)
    {
        regionName = ".";
    }

    return fileName(commsDir_/regionName/this->patch().name());
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::writeGeometry() const
{
    int tag = Pstream::msgType() + 1;

    const label procI = Pstream::myProcNo();
    const polyPatch& p = this->patch().patch();
    const polyMesh& mesh = p.boundaryMesh().mesh();

    labelList pointToGlobal;
    labelList uniquePointIDs;
    (void)mesh.globalData().mergePoints
    (
        p.meshPoints(),
        p.meshPointMap(),
        pointToGlobal,
        uniquePointIDs
    );

    List<pointField> allPoints(Pstream::nProcs());
    allPoints[procI] = pointField(mesh.points(), uniquePointIDs);
    Pstream::gatherList(allPoints, tag);

    List<faceList> allFaces(Pstream::nProcs());
    faceList& patchFaces = allFaces[procI];
    patchFaces = p.localFaces();
    forAll(patchFaces, faceI)
    {
        inplaceRenumber(pointToGlobal, patchFaces[faceI]);
    }

    Pstream::gatherList(allFaces, tag);

    if (Pstream::master())
    {
        OFstream osPoints(baseDir()/"patchPoints");
        if (log_)
        {
            Info<< "writing patch points to: " << osPoints.name() << endl;
        }

        osPoints<<
            ListListOps::combine<pointField>(allPoints, accessOp<pointField>());

        OFstream osFaces(baseDir()/"patchFaces");
        if (log_)
        {
            Info<< "writing patch faces to: " << osFaces.name() << endl;
        }

        osFaces<<
            ListListOps::combine<faceList>(allFaces, accessOp<faceList>());
    }
}


template<class Type>
Foam::fileName Foam::externalCoupledMixedFvPatchField<Type>::lockFile() const
{
    return fileName(baseDir()/(lockName + ".lock"));
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::createLockFile() const
{
    if (!Pstream::master())
    {
        return;
    }

    if (log_)
    {
        Info<< type() << ": creating lock file" << endl;
    }

    OFstream os(lockFile());
    os  << "waiting";
    os.flush();
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::removeLockFile() const
{
    if (!Pstream::master())
    {
        return;
    }

    if (log_)
    {
        Info<< type() << ": removing lock file" << endl;
    }

    rm(lockFile());
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::writeAndWait
(
    const fileName& transferFile
) const
{
    if (log_)
    {
        Info<< type() << ": writing data to " << transferFile << endl;
    }

    if (Pstream::parRun())
    {
        int tag = Pstream::msgType() + 1;

        List<Field<Type> > values(Pstream::nProcs());
        values[Pstream::myProcNo()].setSize(this->refValue().size());
        values[Pstream::myProcNo()] = this->refValue();
        Pstream::gatherList(values, tag);

        List<Field<Type> > grads(Pstream::nProcs());
        grads[Pstream::myProcNo()].setSize(this->refGrad().size());
        grads[Pstream::myProcNo()] = this->refGrad();
        Pstream::gatherList(grads, tag);

        List<scalarField> fracs(Pstream::nProcs());
        fracs[Pstream::myProcNo()].setSize(this->valueFraction().size());
        fracs[Pstream::myProcNo()] = this->valueFraction();
        Pstream::gatherList(fracs, tag);

        if (Pstream::master())
        {
            OFstream os(transferFile);

            forAll(values, procI)
            {
                const Field<Type>& v = values[procI];
                const Field<Type>& g = grads[procI];
                const scalarField& f = fracs[procI];

                forAll(v, faceI)
                {
                    os  << v[faceI] << token::SPACE
                        << g[faceI] << token::SPACE
                        << f[faceI] << nl;
                }
            }

            os.flush();
        }
    }
    else
    {
        OFstream os(transferFile);

        forAll(this->patch(), faceI)
        {
            os  << this->refValue()[faceI] << token::SPACE
                << this->refGrad()[faceI] << token::SPACE
                << this->valueFraction()[faceI] << nl;
        }

        os.flush();
    }

    // remove lock file, signalling external source to execute
    removeLockFile();


    if (log_)
    {
        Info<< type() << ": beginning wait for lock file " << lockFile()
            << endl;
    }

    bool found = false;
    label totalTime = 0;

    while (!found)
    {
        sleep(waitInterval_);
        totalTime += waitInterval_;

        if (log_)
        {
            Info<< type() << ": wait time = " << totalTime << endl;
        }

        if (totalTime > timeOut_)
        {
            FatalErrorIn
            (
                "void Foam::externalCoupledMixedFvPatchField<Type>::"
                "writeAndWait(const fileName&) const"
            )
                << "Wait time exceeded time out time of " << timeOut_
                << " s" << abort(FatalError);
        }

        IFstream is(lockFile());

        if (is.good())
        {
            if (log_)
            {
                Info<< type() << ": found lock file " << lockFile() << endl;
            }

            found = true;
        }
    }


    if (Pstream::master())
    {
        // remove old data file from OpenFOAM
        rm(transferFile);
    }
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::initialiseRead
(
    IFstream& is
) const
{
    if (!is.good())
    {
        FatalErrorIn
        (
            "void Foam::externalCoupledMixedFvPatchField<Type>::"
            "initialiseRead()"
        )
            << "Unable to open data transfer file " << is.name()
            << " for patch " << this->patch().name()
            << exit(FatalError);
    }

    if (Pstream::parRun())
    {
        // fast-forward to relevant point in file
        globalIndex gi(this->patch().size());

        if (this->patch().size())
        {
            string line;
            const label offset = gi.offset(Pstream::myProcNo());
            for (label i = 0; i < offset; i++)
            {
                if (is.good())
                {
                    is.getLine(line);
                }
                else
                {
                    FatalErrorIn
                    (
                        "void Foam::externalCoupledMixedFvPatchField<Type>::"
                        "initialiseRead()"
                    )
                        << "Unable to distribute parallel data for file "
                        << is.name() << " for patch " << this->patch().name()
                        << exit(FatalError);
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::externalCoupledMixedFvPatchField<Type>::externalCoupledMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    commsDir_("unknown-commsDir"),
    waitInterval_(0),
    timeOut_(0),
    calcFrequency_(0),
    log_(false)
{
    this->refValue() = pTraits<Type>::zero;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::externalCoupledMixedFvPatchField<Type>::externalCoupledMixedFvPatchField
(
    const externalCoupledMixedFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    commsDir_(ptf.commsDir_),
    fName_(ptf.fName_),
    waitInterval_(ptf.waitInterval_),
    timeOut_(ptf.timeOut_),
    calcFrequency_(ptf.calcFrequency_),
    log_(ptf.log_)
{}


template<class Type>
Foam::externalCoupledMixedFvPatchField<Type>::externalCoupledMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
    commsDir_(dict.lookup("commsDir")),
    fName_(dict.lookup("fileName")),
    waitInterval_(dict.lookupOrDefault("waitInterval", 1)),
    timeOut_(dict.lookupOrDefault("timeOut", 100*waitInterval_)),
    calcFrequency_(dict.lookupOrDefault("calcFrequency", 1)),
    log_(dict.lookupOrDefault("log", false))
{
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<Type>::operator=(this->patchInternalField());
    }

    if (Pstream::master())
    {
        commsDir_.expand();
        mkDir(baseDir());
        createLockFile();
    }

    // initialise as a fixed value
    this->refValue() = *this;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 1.0;

    writeGeometry();
}


template<class Type>
Foam::externalCoupledMixedFvPatchField<Type>::externalCoupledMixedFvPatchField
(
    const externalCoupledMixedFvPatchField& ecmpf
)
:
    mixedFvPatchField<Type>(ecmpf),
    commsDir_(ecmpf.commsDir_),
    fName_(ecmpf.fName_),
    waitInterval_(ecmpf.waitInterval_),
    timeOut_(ecmpf.timeOut_),
    calcFrequency_(ecmpf.calcFrequency_),
    log_(ecmpf.log_)
{}


template<class Type>
Foam::externalCoupledMixedFvPatchField<Type>::externalCoupledMixedFvPatchField
(
    const externalCoupledMixedFvPatchField& ecmpf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ecmpf, iF),
    commsDir_(ecmpf.commsDir_),
    fName_(ecmpf.fName_),
    waitInterval_(ecmpf.waitInterval_),
    timeOut_(ecmpf.timeOut_),
    calcFrequency_(ecmpf.calcFrequency_),
    log_(ecmpf.log_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (this->db().time().timeIndex() % calcFrequency_ == 0)
    {
        fileName transferFile(baseDir()/fName_);

        // write data for external source and wait for response
        writeAndWait(transferFile + ".out");

        // read data passed back from external source
        IFstream is(transferFile + ".in");

        // pre-process the input transfer file
        initialiseRead(is);

        // read data from file
        forAll(this->patch(), faceI)
        {
            if (is.good())
            {
                is  >> this->refValue()[faceI]
                    >> this->refGrad()[faceI]
                    >> this->valueFraction()[faceI];
            }
            else
            {
                FatalErrorIn
                (
                    "void Foam::externalCoupledMixedFvPatchField<Type>::"
                    "updateCoeffs()"
                )
                    << "Insufficient data for patch " << this->patch().name()
                    << " in file " << is.name() << exit(FatalError);
            }
        }

        // create lock file for external source
        createLockFile();
    }

    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::write(Ostream& os) const
{
    mixedFvPatchField<Type>::write(os);

    os.writeKeyword("commsDir") << commsDir_ << token::END_STATEMENT << nl;
    os.writeKeyword("fileName") << fName_ << token::END_STATEMENT << nl;
    os.writeKeyword("waitInterval") << waitInterval_ << token::END_STATEMENT
        << nl;
    os.writeKeyword("timeOut") << timeOut_ << token::END_STATEMENT << nl;
    os.writeKeyword("calcFrequency") << calcFrequency_ << token::END_STATEMENT
        << nl;
    os.writeKeyword("log") << log_ << token::END_STATEMENT << nl;

    this->writeEntry("value", os);
}


// ************************************************************************* //
