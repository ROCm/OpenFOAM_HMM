/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "fileOperation.H"
#include "uncollatedFileOperation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::refPtr<Foam::fileOperation> Foam::fileOperation::fileHandlerPtr_;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::fileOperation& Foam::fileOperation::fileHandler()
{
    if (!fileOperation::fileHandlerPtr_)
    {
        word handlerType(Foam::getEnv("FOAM_FILEHANDLER"));

        if (handlerType.empty())
        {
            handlerType = defaultFileHandler;
        }

        fileOperation::fileHandlerPtr_ = fileOperation::New(handlerType, true);
    }

    return *fileOperation::fileHandlerPtr_;
}


Foam::refPtr<Foam::fileOperation>
Foam::fileOperation::fileHandler(std::nullptr_t)
{
    return refPtr<fileOperation>(std::move(fileHandlerPtr_));
}


Foam::refPtr<Foam::fileOperation>
Foam::fileOperation::fileHandler(refPtr<fileOperation>& newHandler)
{
    // - do nothing if newHandler is empty. Does not delete current
    // - do nothing if newHandler is identical to current handler

    // Change ownership as atomic operations

    // If newHandler and current handler are actually identical, we
    // have a bit problem somewhere else since this means that the pointer
    // is managed is done in two places!
    // Should flag as a FatalError (in the future), but there may still be
    // some place where we would like to fake shared pointers?

    // TBD: add a flush() operation on the old handler first,
    // instead of waiting for it to be run on destruction?

    refPtr<fileOperation> old;

    if
    (
        newHandler.get() != nullptr
     && newHandler.get() != fileOperation::fileHandlerPtr_.get()
    )
    {
        old.swap(newHandler);
        old.swap(fileOperation::fileHandlerPtr_);
    }

    return old;
}


Foam::refPtr<Foam::fileOperation>
Foam::fileOperation::fileHandler(autoPtr<fileOperation>&& newHandler)
{
    // Same logic as refPtr version

    refPtr<fileOperation> old;

    if
    (
        newHandler.get() != nullptr
     && newHandler.get() != fileOperation::fileHandlerPtr_.get()
    )
    {
        old.reset(newHandler.release());
        old.swap(fileOperation::fileHandlerPtr_);
    }

    return old;
}


Foam::refPtr<Foam::fileOperation>
Foam::fileOperation::fileHandler(refPtr<fileOperation>&& newHandler)
{
    return fileOperation::fileHandler(newHandler);
}


Foam::autoPtr<Foam::fileOperation> Foam::fileOperation::NewUncollated()
{
    return autoPtr<fileOperation>
    (
        new fileOperations::uncollatedFileOperation(false)
    );
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fileOperation>
Foam::fileOperation::New
(
    const word& handlerType,
    bool verbose
)
{
    if (handlerType.empty())
    {
        if (fileOperation::defaultFileHandler.empty())
        {
            FatalErrorInFunction
                << "Default file-handler name is undefined" << nl
                << abort(FatalError);
        }

        return fileOperation::New(fileOperation::defaultFileHandler, verbose);
    }

    DebugInFunction
        << "Constructing fileHandler" << endl;

    auto* ctorPtr = wordConstructorTable(handlerType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "fileHandler",
            handlerType,
            *wordConstructorTablePtr_
        ) << abort(FatalError);
    }

    return autoPtr<fileOperation>(ctorPtr(verbose));
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fileOperation>
Foam::fileHandler(autoPtr<fileOperation>&& newHandler)
{
    refPtr<fileOperation> oldHandler
    (
        fileOperation::fileHandler(std::move(newHandler))
    );

    autoPtr<fileOperation> old;

    // Can return as autoPtr if handler was also a pointer (not a reference)
    if (oldHandler.is_pointer())
    {
        old.reset(oldHandler.release());
    }

    return old;
}


// ************************************************************************* //
