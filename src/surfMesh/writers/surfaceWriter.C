/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "surfaceWriter.H"
#include "proxySurfaceWriter.H"
#include "MeshedSurfaceProxy.H"

#include "Time.H"
#include "globalIndex.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceWriter, 0);
    defineRunTimeSelectionTable(surfaceWriter, word);
    defineRunTimeSelectionTable(surfaceWriter, wordDict);
}

Foam::scalar Foam::surfaceWriter::defaultMergeDim = 1e-8;

const Foam::meshedSurf::emptySurface Foam::surfaceWriter::emptySurface_;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::surfaceWriter::supportedType(const word& writeType)
{
    return
    (
        wordConstructorTablePtr_->found(writeType)
     || wordDictConstructorTablePtr_->found(writeType)
     || MeshedSurfaceProxy<face>::canWriteType(writeType)
    );
}


Foam::autoPtr<Foam::surfaceWriter>
Foam::surfaceWriter::New(const word& writeType)
{
    // Constructors without dictionary options
    auto* ctorPtr = wordConstructorTable(writeType);

    if (!ctorPtr)
    {
        if (MeshedSurfaceProxy<face>::canWriteType(writeType))
        {
            // Generally unknown, but handle via 'proxy' handler
            return autoPtr<surfaceWriter>
            (
                new surfaceWriters::proxyWriter(writeType)
            );
        }

        FatalErrorInFunction
            << "Unknown write type \"" << writeType << "\"\n\n"
            << "Valid write types : "
            << flatOutput(wordConstructorTablePtr_->sortedToc()) << nl
            << "Valid proxy types : "
            << MeshedSurfaceProxy<face>::writeTypes() << endl
            << exit(FatalError);
    }

    return autoPtr<surfaceWriter>(ctorPtr());
}


Foam::autoPtr<Foam::surfaceWriter>
Foam::surfaceWriter::New
(
    const word& writeType,
    const dictionary& writeOpts
)
{
    // Constructors with dictionary options
    {
        auto* ctorPtr = wordDictConstructorTable(writeType);

        if (ctorPtr)
        {
            return autoPtr<surfaceWriter>(ctorPtr(writeOpts));
        }
    }


    // Constructors without dictionary options
    auto* ctorPtr = wordConstructorTable(writeType);

    if (!ctorPtr)
    {
        if (MeshedSurfaceProxy<face>::canWriteType(writeType))
        {
            // Generally unknown, but handle via 'proxy' handler
            return autoPtr<surfaceWriter>
            (
                new surfaceWriters::proxyWriter(writeType, writeOpts)
            );
        }

        FatalErrorInFunction
            << "Unknown write type \"" << writeType << "\"\n\n"
            << "Valid write types : "
            << wordConstructorTablePtr_->sortedToc() << nl
            << "Valid proxy types : "
            << MeshedSurfaceProxy<face>::writeTypes() << endl
            << exit(FatalError);
    }

    return autoPtr<surfaceWriter>(ctorPtr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriter::surfaceWriter()
:
    surf_(std::cref<meshedSurf>(emptySurface_)),
    surfComp_(),
    useComponents_(false),
    upToDate_(false),
    wroteGeom_(false),
    parallel_(true),
    useTimeDir_(false),
    isPointData_(false),
    verbose_(false),
    nFields_(0),
    mergeDim_(defaultMergeDim),
    merged_(),
    currTime_(),
    outputPath_()
{
    surfaceWriter::close();
}


Foam::surfaceWriter::surfaceWriter(const dictionary& options)
:
    surfaceWriter()
{
    options.readIfPresent("verbose", verbose_);
}


Foam::surfaceWriter::surfaceWriter
(
    const meshedSurf& surf,
    bool parallel,
    const dictionary& options
)
:
    surfaceWriter(options)
{
    setSurface(surf, parallel);
}


Foam::surfaceWriter::surfaceWriter
(
    const pointField& points,
    const faceList& faces,
    bool parallel,
    const dictionary& options
)
:
    surfaceWriter(options)
{
    setSurface(points, faces, parallel);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceWriter::~surfaceWriter()
{
    close();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfaceWriter::setTime(const instant& inst)
{
    currTime_ = inst;
}


void Foam::surfaceWriter::setTime(scalar timeValue)
{
    currTime_ = instant(timeValue);
}


void Foam::surfaceWriter::setTime(scalar timeValue, const word& timeName)
{
    currTime_.value() = timeValue;
    currTime_.name() = timeName;
}


void Foam::surfaceWriter::unsetTime()
{
    currTime_.value() = 0;
    currTime_.name().clear();
}


void Foam::surfaceWriter::beginTime(const Time& t)
{
    setTime(t.value(), t.timeName());
}


void Foam::surfaceWriter::beginTime(const instant& inst)
{
    setTime(inst);
}


void Foam::surfaceWriter::endTime()
{
    unsetTime();
}


void Foam::surfaceWriter::open(const fileName& outputPath)
{
    outputPath_ = outputPath;
    wroteGeom_ = false;
}


void Foam::surfaceWriter::open
(
    const meshedSurf& surf,
    const fileName& outputPath,
    bool parallel
)
{
    close();
    setSurface(surf, parallel);
    open(outputPath);
}


void Foam::surfaceWriter::open
(
    const pointField& points,
    const faceList& faces,
    const fileName& outputPath,
    bool parallel
)
{
    close();
    setSurface(points, faces, parallel);
    open(outputPath);
}


void Foam::surfaceWriter::open
(
    const meshedSurf& surf,
    const fileName& outputPath
)
{
    close();
    setSurface(surf, parallel_);
    open(outputPath);
}


void Foam::surfaceWriter::open
(
    const pointField& points,
    const faceList& faces,
    const fileName& outputPath
)
{
    close();
    setSurface(points, faces, parallel_);
    open(outputPath);
}


void Foam::surfaceWriter::close()
{
    outputPath_.clear();
    wroteGeom_ = false;
}


void Foam::surfaceWriter::clear()
{
    close();
    expire();
    useComponents_ = false;
    surf_ = std::cref<meshedSurf>(emptySurface_);
    surfComp_.clear();
}


void Foam::surfaceWriter::setSurface
(
    const meshedSurf& surf,
    bool parallel
)
{
    expire();
    useComponents_ = false;
    surf_ = std::cref<meshedSurf>(surf);
    surfComp_.clear();
    parallel_ = (parallel && Pstream::parRun());
}


void Foam::surfaceWriter::setSurface
(
    const pointField& points,
    const faceList& faces,
    bool parallel
)
{
    expire();
    useComponents_ = true;
    surf_ = std::cref<meshedSurf>(emptySurface_);
    surfComp_.reset(points, faces);
    parallel_ = (parallel && Pstream::parRun());
}


void Foam::surfaceWriter::setSurface
(
    const meshedSurf& surf
)
{
    setSurface(surf, parallel_);
}


void Foam::surfaceWriter::setSurface
(
    const pointField& points,
    const faceList& faces
)
{
    setSurface(points, faces, parallel_);
}


bool Foam::surfaceWriter::needsUpdate() const
{
    return !upToDate_;
}


bool Foam::surfaceWriter::wroteData() const
{
    return wroteGeom_;
}


bool Foam::surfaceWriter::expire()
{
    const bool changed = upToDate_;

    upToDate_ = false;
    wroteGeom_ = false;
    merged_.clear();

    // Field count (nFields_) is a different type of accounting
    // and is unaffected by geometry changes

    return changed;
}


bool Foam::surfaceWriter::hasSurface() const
{
    return (useComponents_ || (&emptySurface_ != &(surf_.get())));
}


bool Foam::surfaceWriter::empty() const
{
    const bool value =
    (
        useComponents_
      ? surfComp_.faces().empty()
      : surf_.get().faces().empty()
    );

    return (parallel_ ? returnReduce(value, andOp<bool>()) : value);
}


Foam::label Foam::surfaceWriter::size() const
{
    const label value =
    (
        useComponents_
      ? surfComp_.faces().size()
      : surf_.get().faces().size()
    );

    return (parallel_ ? returnReduce(value, sumOp<label>()) : value);
}


bool Foam::surfaceWriter::checkOpen() const
{
    if (outputPath_.empty())
    {
        FatalErrorInFunction
            << type() << " : Attempted to write without a path" << nl
            << exit(FatalError);
    }

    return !outputPath_.empty();
}


bool Foam::surfaceWriter::merge() const
{
    bool changed = false;

    if (parallel_ && Pstream::parRun() && !upToDate_)
    {
        if (useComponents_)
        {
            changed = merged_.merge(surfComp_, mergeDim_);
        }
        else
        {
            changed = merged_.merge(surf_.get(), mergeDim_);
        }
    }
    upToDate_ = true;

    if (changed)
    {
        wroteGeom_ = false;
    }

    return changed;
}


const Foam::meshedSurf& Foam::surfaceWriter::surface() const
{
    merge();

    if (parallel_ && Pstream::parRun())
    {
        return merged_;
    }

    if (useComponents_)
    {
        return surfComp_;
    }
    else
    {
        return surf_.get();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::surfaceWriter::mergeFieldTemplate
(
    const Field<Type>& fld
) const
{
    if (parallel_ && Pstream::parRun())
    {
        // Ensure geometry is also merged
        merge();

        // Gather all values
        auto tfield = tmp<Field<Type>>::New();
        auto& allFld = tfield.ref();

        globalIndex::gatherOp(fld, allFld);

        // Renumber (point data) to correspond to merged points
        if
        (
            Pstream::master()
         && this->isPointData()
         && merged_.pointsMap().size()
        )
        {
            inplaceReorder(merged_.pointsMap(), allFld);
            allFld.resize(merged_.points().size());
        }

        return tfield;
    }

    // Mark that any geometry changes have been taken care of
    upToDate_ = true;

    return fld;
}


#define defineSurfaceWriterMergeMethod(ThisClass, Type)                        \
    Foam::tmp<Foam::Field<Type>>                                               \
    ThisClass::mergeField(const Field<Type>& fld) const                        \
    {                                                                          \
        return mergeFieldTemplate(fld);                                        \
    }

defineSurfaceWriterMergeMethod(Foam::surfaceWriter, Foam::label);
defineSurfaceWriterMergeMethod(Foam::surfaceWriter, Foam::scalar);
defineSurfaceWriterMergeMethod(Foam::surfaceWriter, Foam::vector);
defineSurfaceWriterMergeMethod(Foam::surfaceWriter, Foam::sphericalTensor);
defineSurfaceWriterMergeMethod(Foam::surfaceWriter, Foam::symmTensor);
defineSurfaceWriterMergeMethod(Foam::surfaceWriter, Foam::tensor)

#undef defineSurfaceWriterMergeMethod


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<surfaceWriter>& ip
)
{
    const surfaceWriter& w = ip.t_;

    os  << "surfaceWriter:"
        << " upToDate: " << w.upToDate_
        << " PointData: " << w.isPointData_
        << " nFields: " << w.nFields_
        << " time: " << w.currTime_
        << " path: " << w.outputPath_ << endl;

    return os;
}


// ************************************************************************* //
