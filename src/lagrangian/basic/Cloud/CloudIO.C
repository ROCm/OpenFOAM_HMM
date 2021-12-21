/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017, 2020 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "Cloud.H"
#include "Time.H"
#include "IOPosition.H"
#include "IOdictionary.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParticleType>
Foam::word Foam::Cloud<ParticleType>::cloudPropertiesName("cloudProperties");


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class ParticleType>
void Foam::Cloud<ParticleType>::readCloudUniformProperties()
{
    IOobject dictObj
    (
        cloudPropertiesName,
        time().timeName(),
        "uniform"/cloud::prefix/name(),
        db(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        false
    );

    if (dictObj.typeHeaderOk<IOdictionary>(true))
    {
        const IOdictionary uniformPropsDict(dictObj);

        // Fall back to positions mode if the entry is not present for
        // backwards compatibility
        geometryType_ =
            cloud::geometryTypeNames.getOrDefault
            (
                "geometry",
                uniformPropsDict,
                cloud::geometryType::POSITIONS
            );

        const word procName("processor" + Foam::name(Pstream::myProcNo()));

        const dictionary* dictptr = uniformPropsDict.findDict(procName);

        if (dictptr)
        {
            dictptr->readEntry("particleCount", ParticleType::particleCount_);
        }
    }
    else
    {
        ParticleType::particleCount_ = 0;
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::writeCloudUniformProperties() const
{
    IOdictionary uniformPropsDict
    (
        IOobject
        (
            cloudPropertiesName,
            time().timeName(),
            "uniform"/cloud::prefix/name(),
            db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    labelList np(Pstream::nProcs(), Zero);
    np[Pstream::myProcNo()] = ParticleType::particleCount_;

    Pstream::listCombineGather(np, maxEqOp<label>());
    Pstream::listCombineScatter(np);

    uniformPropsDict.add
    (
        "geometry",
        cloud::geometryTypeNames[geometryType_]
    );

    forAll(np, i)
    {
        word procName("processor" + Foam::name(i));
        uniformPropsDict.add(procName, dictionary());
        uniformPropsDict.subDict(procName).add("particleCount", np[i]);
    }

    uniformPropsDict.writeObject
    (
        IOstreamOption(IOstream::ASCII, time().writeCompression()),
        true
    );
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::initCloud(const bool checkClass)
{
    readCloudUniformProperties();

    IOPosition<Cloud<ParticleType>> ioP(*this, geometryType_);

    const bool valid = ioP.headerOk();
    Istream& is = ioP.readStream(checkClass ? typeName : "", valid);
    if (valid)
    {
        ioP.readData(is, *this);
        ioP.close();
    }

    if (!valid && debug)
    {
        Pout<< "Cannot read particle positions file:" << nl
            << "    " << ioP.objectPath() << nl
            << "Assuming the initial cloud contains 0 particles." << endl;
    }

    // Always operate in coordinates mode after reading
    geometryType_ = cloud::geometryType::COORDINATES;

    // Ask for the tetBasePtIs to trigger all processors to build
    // them, otherwise, if some processors have no particles then
    // there is a comms mismatch.
    polyMesh_.tetBasePtIs();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Cloud<ParticleType>::Cloud
(
    const polyMesh& pMesh,
    const word& cloudName,
    const bool checkClass
)
:
    cloud(pMesh, cloudName),
    polyMesh_(pMesh),
    labels_(),
    cellWallFacesPtr_(),
    geometryType_(cloud::geometryType::COORDINATES)
{
    checkPatches();

    polyMesh_.tetBasePtIs();
    polyMesh_.oldCellCentres();

    initCloud(checkClass);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
Foam::IOobject Foam::Cloud<ParticleType>::fieldIOobject
(
    const word& fieldName,
    const IOobject::readOption r
) const
{
    return IOobject
    (
        fieldName,
        time().timeName(),
        *this,
        r,
        IOobject::NO_WRITE,
        false
    );
}


template<class ParticleType>
template<class DataType>
void Foam::Cloud<ParticleType>::checkFieldIOobject
(
    const Cloud<ParticleType>& c,
    const IOField<DataType>& data
) const
{
    if (data.size() != c.size())
    {
        FatalErrorInFunction
            << "Size of " << data.name()
            << " field " << data.size()
            << " does not match the number of particles " << c.size()
            << abort(FatalError);
    }
}


template<class ParticleType>
template<class DataType>
void Foam::Cloud<ParticleType>::checkFieldFieldIOobject
(
    const Cloud<ParticleType>& c,
    const CompactIOField<Field<DataType>, DataType>& data
) const
{
    if (data.size() != c.size())
    {
        FatalErrorInFunction
            << "Size of " << data.name()
            << " field " << data.size()
            << " does not match the number of particles " << c.size()
            << abort(FatalError);
    }
}


template<class ParticleType>
template<class Type>
bool Foam::Cloud<ParticleType>::readStoreFile
(
    const IOobject& io,
    const IOobject& ioNew
) const
{
    if (io.headerClassName() == IOField<Type>::typeName)
    {
        IOField<Type> fld(io);
        auto* fldNewPtr = new IOField<Type>(ioNew, std::move(fld));
        return fldNewPtr->store();
    }

    return false;
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::readFromFiles
(
    objectRegistry& obr,
    const wordRes& selectFields
) const
{
    IOobjectList cloudObjects
    (
        *this,
        time().timeName(),
        "",
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    forAllIters(cloudObjects, iter)
    {
        if (selectFields.size() && !selectFields.match(iter()->name()))
        {
            continue;
        }

        IOobject ioNew
        (
            iter()->name(),
            time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        );

        auto& object = *iter();

        const bool stored
        (
            readStoreFile<label>(object, ioNew)
         || readStoreFile<scalar>(object, ioNew)
         || readStoreFile<vector>(object, ioNew)
         || readStoreFile<sphericalTensor>(object, ioNew)
         || readStoreFile<symmTensor>(object, ioNew)
         || readStoreFile<tensor>(object, ioNew)
        );

        if (!stored)
        {
            DebugInfo
                << "Unhandled field type " << iter()->headerClassName()
                << endl;
        }
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::writeFields() const
{
    ParticleType::writeFields(*this);
}


template<class ParticleType>
bool Foam::Cloud<ParticleType>::writeObject
(
    IOstreamOption streamOpt,
    const bool
) const
{
    writeCloudUniformProperties();

    writeFields();
    return cloud::writeObject(streamOpt, this->size());
}


// * * * * * * * * * * * * * * * Ostream Operators * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Ostream& Foam::operator<<(Ostream& os, const Cloud<ParticleType>& c)
{
    c.writeData(os);

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
