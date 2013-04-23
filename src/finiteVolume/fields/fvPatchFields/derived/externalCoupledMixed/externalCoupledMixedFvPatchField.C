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
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "IFstream.H"
#include "globalIndex.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
Foam::word Foam::externalCoupledMixedFvPatchField<Type>::lockName = "OpenFOAM";

template<class Type>
Foam::string
Foam::externalCoupledMixedFvPatchField<Type>::patchKey = "# Patch: ";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::externalCoupledMixedFvPatchField<Type>::baseDir
(
    const word& patchName
) const
{
    word regionName(this->dimensionedInternalField().mesh().name());
    if (regionName == polyMesh::defaultRegion)
    {
        regionName = ".";
    }

    fileName result(commsDir_/regionName);
    result.clean();

    if (collate_)
    {
        return result;
    }
    else
    {
        if (patchName == word::null)
        {
            return fileName(result/this->patch().name());
        }
        else
        {
            return fileName(result/patchName);
        }
    }
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::setMaster()
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;

    const volFieldType& cvf =
        static_cast<const volFieldType&>(this->dimensionedInternalField());

    volFieldType& vf = const_cast<volFieldType&>(cvf);

    typename volFieldType::GeometricBoundaryField& bf = vf.boundaryField();

    if (collate_)
    {
        bool found = false;
        forAll(bf, patchI)
        {
            if (isA<externalCoupledMixedFvPatchField<Type> >(bf[patchI]))
            {
                externalCoupledMixedFvPatchField<Type>& pf =
                    refCast<externalCoupledMixedFvPatchField<Type> >
                    (
                        bf[patchI]
                    );

                // only attempt to change master flags of BCs that have not
                // been set (or at least only the master)
                if (pf.master())
                {
                    if (!found)
                    {
                        pf.master() = true;
                        found = true;
                    }
                    else
                    {
                        pf.master() = false;
                    }
                }
            }
        }
    }
    else
    {
        // check that collated flag is not set on any other patches
        forAll(bf, patchI)
        {
            if (isA<externalCoupledMixedFvPatchField<Type> >(bf[patchI]))
            {
                const externalCoupledMixedFvPatchField<Type>& pf =
                    refCast<const externalCoupledMixedFvPatchField<Type> >
                    (
                        bf[patchI]
                    );

                if (pf.collate())
                {
                    FatalErrorIn
                    (
                        "void Foam::externalCoupledMixedFvPatchField<Type>::"
                        "setMaster()"
                    )   << "All " << type() << " patches should either use "
                        << "collate = true OR false, but not a mix of both"
                        << exit(FatalError);
                }
            }
        }

        master_ = true;
    }
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::writeGeometry
(
    OFstream& osPoints,
    OFstream& osFaces
) const
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
        pointField pts
        (
            ListListOps::combine<pointField>(allPoints, accessOp<pointField>())
        );

        // write points
        osPoints << patchKey.c_str() << this->patch().name() << pts << endl;

        faceList fcs
        (
            ListListOps::combine<faceList>(allFaces, accessOp<faceList>())
        );

        // write faces
        osFaces<< patchKey.c_str() << this->patch().name() << fcs << endl;
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
    if (!master_ || !Pstream::master())
    {
        return;
    }

    const fileName fName(lockFile());
    IFstream is(fName);

    // only create lock file if it doesn't already exist
    if (!is.good())
    {
        if (log_)
        {
            Info<< type() << ": creating lock file" << endl;
        }

        OFstream os(fName);
        os  << "lock file";
        os.flush();
    }
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::removeLockFile() const
{
    if (!master_ || !Pstream::master())
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
void Foam::externalCoupledMixedFvPatchField<Type>::startWait() const
{
    if (collate_)
    {
        // only wait on master patch

        typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;

        const volFieldType& cvf =
            static_cast<const volFieldType&>(this->dimensionedInternalField());

        const typename volFieldType::GeometricBoundaryField& bf =
            cvf.boundaryField();

        forAll(bf, patchI)
        {
            if (isA<externalCoupledMixedFvPatchField<Type> >(bf[patchI]))
            {
                const externalCoupledMixedFvPatchField<Type>& pf =
                    refCast<const externalCoupledMixedFvPatchField<Type> >
                    (
                        bf[patchI]
                    );

                if (pf.master())
                {
                    pf.wait();
                    break;
                }
            }
        }
    }
    else
    {
        wait();
    }
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::wait() const
{
    const fileName fName(lockFile());
    bool found = false;
    label totalTime = 0;

    if (log_)
    {
        Info<< type() << ": beginning wait for lock file " << fName << endl;
    }

    while (!found)
    {
        if (totalTime > timeOut_)
        {
            FatalErrorIn
            (
                "void "
                "Foam::externalCoupledMixedFvPatchField<Type>::wait() const"
            )
                << "Wait time exceeded time out time of " << timeOut_
                << " s" << abort(FatalError);
        }

        IFstream is(fName);

        if (is.good())
        {
            if (log_)
            {
                Info<< type() << ": found lock file " << fName << endl;
            }

            found = true;
            break;
        }

        sleep(waitInterval_);
        totalTime += waitInterval_;

        if (log_)
        {
            Info<< type() << ": wait time = " << totalTime << endl;
        }
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
            << "Unable to open data transfer file " << is.name().caseName()
            << " for patch " << this->patch().name()
            << exit(FatalError);
    }

    string line;

    // scan forward to the line that starts '# Patch: <myPatchName>'
    const string searchStr(patchKey + this->patch().name());

    bool scan = true;
    while (is.good() && scan)
    {
        is.getLine(line);

        if (line.rfind(searchStr) != std::string::npos)
        {
            scan = false;
        }
    }

    if (scan)
    {
        FatalErrorIn
        (
            "void Foam::externalCoupledMixedFvPatchField<Type>::"
            "initialiseRead"
            "("
                "IFstream&"
            ") const"
        )
            << "Unable to find data starting with " << searchStr
            << " in file" << nl
            << "    " << is.name().caseName() << abort(FatalError);
    }

    if (Pstream::parRun())
    {
        // fast-forward to relevant point in file
        globalIndex gi(this->patch().size());

        if (this->patch().size())
        {
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
                        "initialiseRead"
                        "("
                            "IFstream&"
                        ") const"
                    )
                        << "Unable to distribute parallel data for file "
                        << is.name().caseName() << " for patch "
                        << this->patch().name() << exit(FatalError);
                }
            }
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::writeData
(
    const fileName& transferFile
) const
{
    if (!master_)
    {
        return;
    }

    OFstream os(transferFile);

    writeHeader(os);

    if (collate_)
    {
        typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;

        const volFieldType& cvf =
            static_cast<const volFieldType&>(this->dimensionedInternalField());

        volFieldType& vf = const_cast<volFieldType&>(cvf);

        typename volFieldType::GeometricBoundaryField& bf = vf.boundaryField();

        forAll(bf, patchI)
        {
            if (isA<externalCoupledMixedFvPatchField<Type> >(bf[patchI]))
            {
                const externalCoupledMixedFvPatchField<Type>& pf =
                    refCast<const externalCoupledMixedFvPatchField<Type> >
                    (
                        bf[patchI]
                    );

                os  << patchKey.c_str() << pf.patch().name() << nl;
                pf.transferData(os);
            }
        }
    }
    else
    {
        os  << patchKey.c_str() << this->patch().name() << nl;
        transferData(os);
    }
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::writeHeader
(
    OFstream& os
) const
{
    os  << "# Values: magSf value snGrad" << endl;
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
    fName_("unknown-fName"),
    collate_(false),
    waitInterval_(0),
    timeOut_(0),
    calcFrequency_(0),
    log_(false),
    master_(false)
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
    collate_(ptf.collate_),
    waitInterval_(ptf.waitInterval_),
    timeOut_(ptf.timeOut_),
    calcFrequency_(ptf.calcFrequency_),
    log_(ptf.log_),
    master_(ptf.master_)
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
    collate_(readBool(dict.lookup("collate"))),
    waitInterval_(dict.lookupOrDefault("waitInterval", 1)),
    timeOut_(dict.lookupOrDefault("timeOut", 100*waitInterval_)),
    calcFrequency_(dict.lookupOrDefault("calcFrequency", 1)),
    log_(dict.lookupOrDefault("log", false)),
    master_(true)
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
    }

    createLockFile();

    // initialise as a fixed value
    this->refValue() = *this;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 1.0;
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
    collate_(ecmpf.collate_),
    waitInterval_(ecmpf.waitInterval_),
    timeOut_(ecmpf.timeOut_),
    calcFrequency_(ecmpf.calcFrequency_),
    log_(ecmpf.log_),
    master_(ecmpf.master_)
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
    collate_(ecmpf.collate_),
    waitInterval_(ecmpf.waitInterval_),
    timeOut_(ecmpf.timeOut_),
    calcFrequency_(ecmpf.calcFrequency_),
    log_(ecmpf.log_),
    master_(ecmpf.master_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::externalCoupledMixedFvPatchField<Type>::
~externalCoupledMixedFvPatchField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    setMaster();

    if (this->db().time().timeIndex() % calcFrequency_ == 0)
    {
        fileName transferFile(baseDir()/fName_);

        // write data for external source
        writeData(transferFile + ".out");

        // remove lock file, signalling external source to execute
        removeLockFile();

        // wait for response
        startWait();

        if (master_ && Pstream::master())
        {
            // remove old data file from OpenFOAM
            rm(transferFile + ".out");
        }

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
                    << " in file " << is.name().caseName() << exit(FatalError);
            }
        }

        // create lock file for external source
        createLockFile();
    }

    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::transferData
(
    OFstream& os
) const
{
    if (log_)
    {
        Info<< type() << ": writing data to " << os.name().caseName() << endl;
    }

    if (Pstream::parRun())
    {
        int tag = Pstream::msgType() + 1;

        List<Field<scalar> > magSfs(Pstream::nProcs());
        magSfs[Pstream::myProcNo()].setSize(this->patch().size());
        magSfs[Pstream::myProcNo()] = this->patch().magSf();
        Pstream::gatherList(magSfs, tag);

        List<Field<Type> > values(Pstream::nProcs());
        values[Pstream::myProcNo()].setSize(this->patch().size());
        values[Pstream::myProcNo()] = this->refValue();
        Pstream::gatherList(values, tag);

        List<Field<Type> > snGrads(Pstream::nProcs());
        snGrads[Pstream::myProcNo()].setSize(this->patch().size());
        snGrads[Pstream::myProcNo()] = this->snGrad();
        Pstream::gatherList(snGrads, tag);

        if (Pstream::master())
        {
            forAll(values, procI)
            {
                const Field<scalar>& magSf = magSfs[procI];
                const Field<Type>& value = values[procI];
                const Field<Type>& snGrad = snGrads[procI];

                forAll(magSf, faceI)
                {
                    os  << magSf[faceI] << token::SPACE
                        << value[faceI] << token::SPACE
                        << snGrad[faceI] << nl;
                }
            }

            os.flush();
        }
    }
    else
    {
        const Field<scalar>& magSf(this->patch().magSf());
        const Field<Type>& value(this->refValue());
        const Field<Type> snGrad(this->snGrad());

        forAll(magSf, faceI)
        {
            os  << magSf[faceI] << token::SPACE
                << value[faceI] << token::SPACE
                << snGrad[faceI] << nl;
        }

        os.flush();
    }
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::writeGeometry() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;

    const volFieldType& cvf =
        static_cast<const volFieldType&>(this->dimensionedInternalField());

    const typename volFieldType::GeometricBoundaryField& bf =
        cvf.boundaryField();

    if (collate_)
    {
        OFstream osPoints(baseDir()/"patchPoints");
        OFstream osFaces(baseDir()/"patchFaces");

        if (log_)
        {
            Info<< "writing collated patch points to: "
                << osPoints.name().caseName() << endl;
            Info<< "writing collated patch faces to: "
                << osFaces.name().caseName() << endl;
        }

        forAll(bf, patchI)
        {
            if (isA<externalCoupledMixedFvPatchField<Type> >(bf[patchI]))
            {
                const externalCoupledMixedFvPatchField<Type>& pf =
                    refCast<const externalCoupledMixedFvPatchField<Type> >
                    (
                        bf[patchI]
                    );

                pf.writeGeometry(osPoints, osFaces);
            }
        }
    }
    else
    {
        forAll(bf, patchI)
        {
            if (isA<externalCoupledMixedFvPatchField<Type> >(bf[patchI]))
            {
                const word& patchName = this->patch().name();

                OFstream osPoints(baseDir(patchName)/"patchPoints");
                OFstream osFaces(baseDir(patchName)/"patchFaces");

                if (log_)
                {
                    Info<< "writing patch " << patchName << " points to: "
                        << osPoints.name().caseName() << endl;
                    Info<< "writing patch " << patchName << " faces to: "
                        << osFaces.name().caseName() << endl;
                }

                const externalCoupledMixedFvPatchField<Type>& pf =
                    refCast<const externalCoupledMixedFvPatchField<Type> >
                    (
                        bf[patchI]
                    );

                pf.writeGeometry(osPoints, osFaces);
            }
        }
    }
}


template<class Type>
void Foam::externalCoupledMixedFvPatchField<Type>::write(Ostream& os) const
{
    mixedFvPatchField<Type>::write(os);

    os.writeKeyword("commsDir") << commsDir_ << token::END_STATEMENT << nl;
    os.writeKeyword("fileName") << fName_ << token::END_STATEMENT << nl;
    os.writeKeyword("collate") << collate_ << token::END_STATEMENT << nl;
    os.writeKeyword("waitInterval") << waitInterval_ << token::END_STATEMENT
        << nl;
    os.writeKeyword("timeOut") << timeOut_ << token::END_STATEMENT << nl;
    os.writeKeyword("calcFrequency") << calcFrequency_ << token::END_STATEMENT
        << nl;
    os.writeKeyword("log") << log_ << token::END_STATEMENT << nl;

    this->writeEntry("value", os);
}


// ************************************************************************* //
