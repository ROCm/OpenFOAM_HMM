/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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

#include "regIOobject.H"
#include "IFstream.H"
#include "Time.H"
#include "Pstream.H"
#include "HashSet.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::regIOobject::read
(
    const bool masterOnly,
    const IOstream::streamFormat format,
    const word& typeName
)
{
    bool ok = true;
    if (Pstream::master() || !masterOnly)
    {
        if (IFstream::debug)
        {
            Pout<< "regIOobject::read() : "
                << "reading object " << name()
                << " from file " << endl;
        }

        // Set flag for e.g. codeStream
        const bool oldGlobal = globalObject();
        globalObject() = masterOnly;
        // If codeStream originates from dictionary which is
        // not IOdictionary we have a problem so use global
        const bool oldFlag = regIOobject::masterOnlyReading;
        regIOobject::masterOnlyReading = masterOnly;

        // Read file
        ok = readData(readStream(typeName));
        close();

        // Restore flags
        globalObject() = oldGlobal;
        regIOobject::masterOnlyReading = oldFlag;
    }

    if (masterOnly && Pstream::parRun())
    {
        // Scatter master data using communication scheme

        const List<Pstream::commsStruct>& comms =
        (
            (Pstream::nProcs() < Pstream::nProcsSimpleSum)
          ? Pstream::linearCommunication()
          : Pstream::treeCommunication()
        );

        // Master reads headerclassname from file. Make sure this gets
        // transfered as well as contents.
        Pstream::scatter
        (
            comms,
            const_cast<word&>(headerClassName()),
            Pstream::msgType(),
            Pstream::worldComm
        );
        Pstream::scatter(comms, note(), Pstream::msgType(), Pstream::worldComm);


        // Get my communication order
        const Pstream::commsStruct& myComm = comms[Pstream::myProcNo()];

        // Reveive from up
        if (myComm.above() != -1)
        {
            if (IFstream::debug)
            {
                Pout<< "regIOobject::read() : "
                    << "reading object " << name()
                    << " from processor " << myComm.above()
                    << endl;
            }

            IPstream fromAbove
            (
                Pstream::scheduled,
                myComm.above(),
                0,
                Pstream::msgType(),
                Pstream::worldComm,
                format
            );
            ok = readData(fromAbove);
        }

        // Send to my downstairs neighbours
        forAll(myComm.below(), belowI)
        {
            OPstream toBelow
            (
                Pstream::scheduled,
                myComm.below()[belowI],
                0,
                Pstream::msgType(),
                Pstream::worldComm,
                format
            );
            bool okWrite = writeData(toBelow);
            ok = ok && okWrite;
        }
    }
    return ok;
}


bool Foam::regIOobject::readHeaderOk
(
    const IOstream::streamFormat format,
    const word& typeName
)
{
    // Everyone check or just master
    bool masterOnly =
        global()
     && (
            regIOobject::fileModificationChecking == timeStampMaster
         || regIOobject::fileModificationChecking == inotifyMaster
        );


    // Check if header is ok for READ_IF_PRESENT
    bool isHeaderOk = false;
    if (readOpt() == IOobject::READ_IF_PRESENT)
    {
        if (masterOnly)
        {
            if (Pstream::master())
            {
                isHeaderOk = headerOk();
            }
            Pstream::scatter(isHeaderOk);
        }
        else
        {
            isHeaderOk = headerOk();
        }
    }

    if
    (
        (
            readOpt() == IOobject::MUST_READ
         || readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || isHeaderOk
    )
    {
        return regIOobject::read(masterOnly, format, typeName);
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Istream& Foam::regIOobject::readStream()
{
    if (IFstream::debug)
    {
        Pout<< "regIOobject::readStream() : "
            << "reading object " << name()
            << " (global " << global() << ")"
            << " from file " << objectPath()
            << endl;
    }

    if (readOpt() == NO_READ)
    {
        FatalErrorInFunction
            << "NO_READ specified for read-constructor of object " << name()
            << " of class " << headerClassName()
            << abort(FatalError);
    }

    // Construct object stream and read header if not already constructed
    if (!isPtr_)
    {
        fileName objPath;
        if (watchIndices_.size())
        {
            // File is being watched. Read exact file that is being watched.
            objPath = time().getFile(watchIndices_.last());
        }
        else
        {
            // Search intelligently for file
            objPath = filePath();

            if (IFstream::debug)
            {
                Pout<< "regIOobject::readStream() : "
                    << "found object " << name()
                    << " (global " << global() << ")"
                    << " in file " << objPath
                    << endl;
            }

            if (!objPath.size())
            {
                FatalIOError
                (
                    "regIOobject::readStream()",
                    __FILE__,
                    __LINE__,
                    objectPath(),
                    0
                )   << "cannot find file"
                    << exit(FatalIOError);
            }
        }

        if (!(isPtr_ = IOobject::objectStream(objPath)))
        {
            FatalIOError
            (
                "regIOobject::readStream()",
                __FILE__,
                __LINE__,
                objPath,
                0
            )   << "cannot open file"
                << exit(FatalIOError);
        }
        else if (!readHeader(*isPtr_))
        {
            FatalIOErrorInFunction(*isPtr_)
                << "problem while reading header for object " << name()
                << exit(FatalIOError);
        }
    }

    return *isPtr_;
}


Foam::Istream& Foam::regIOobject::readStream(const word& expectName)
{
    if (IFstream::debug)
    {
        Pout<< "regIOobject::readStream(const word&) : "
            << "reading object " << name()
            << " of type " << type()
            << " from file " << objectPath()
            << endl;
    }

    // Construct IFstream if not already constructed
    if (!isPtr_)
    {
        readStream();

        // Check the className of the regIOobject
        // dictionary is an allowable name in case the actual class
        // instantiated is a dictionary
        if
        (
            expectName.size()
         && headerClassName() != expectName
         && headerClassName() != "dictionary"
        )
        {
            FatalIOErrorInFunction(*isPtr_)
                << "unexpected class name " << headerClassName()
                << " expected " << expectName << endl
                << "    while reading object " << name()
                << exit(FatalIOError);
        }
    }

    return *isPtr_;
}


void Foam::regIOobject::close()
{
    if (IFstream::debug)
    {
        Pout<< "regIOobject::close() : "
            << "finished reading " << filePath()
            << endl;
    }

    if (isPtr_)
    {
        delete isPtr_;
        isPtr_ = NULL;
    }
}


bool Foam::regIOobject::readData(Istream&)
{
    return false;
}


bool Foam::regIOobject::read()
{
    // Note: cannot do anything in readStream itself since this is used by
    // e.g. GeometricField.


    // Save old watchIndices and clear (so the list of included files can
    // change)
    fileNameList oldWatchFiles;
    if (watchIndices_.size())
    {
        oldWatchFiles.setSize(watchIndices_.size());
        forAll(watchIndices_, i)
        {
            oldWatchFiles[i] = time().getFile(watchIndices_[i]);
        }
        forAllReverse(watchIndices_, i)
        {
            time().removeWatch(watchIndices_[i]);
        }
        watchIndices_.clear();
    }


    // Read
    bool masterOnly =
        global()
     && (
            regIOobject::fileModificationChecking == timeStampMaster
         || regIOobject::fileModificationChecking == inotifyMaster
        );

    bool ok = read(masterOnly, IOstream::BINARY, type());

    if (oldWatchFiles.size())
    {
        // Re-watch master file
        addWatch();
    }

    return ok;
}


bool Foam::regIOobject::modified() const
{
    forAllReverse(watchIndices_, i)
    {
        if (time().getState(watchIndices_[i]) != fileMonitor::UNMODIFIED)
        {
            return true;
        }
    }

    return false;
}


bool Foam::regIOobject::readIfModified()
{
    // Get index of modified file so we can give nice message. Could instead
    // just call above modified()
    label modified = -1;
    forAllReverse(watchIndices_, i)
    {
        if (time().getState(watchIndices_[i]) != fileMonitor::UNMODIFIED)
        {
            modified = watchIndices_[i];
            break;
        }
    }

    if (modified != -1)
    {
        const fileName& fName = time().getFile(watchIndices_.last());

        if (modified == watchIndices_.last())
        {
            Info<< "regIOobject::readIfModified() : " << nl
                << "    Re-reading object " << name()
                << " from file " << fName << endl;
        }
        else
        {
            Info<< "regIOobject::readIfModified() : " << nl
                << "    Re-reading object " << name()
                << " from file " << fName
                << " because of modified file " << time().getFile(modified)
                << endl;
        }

        return read();
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
