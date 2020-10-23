/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

Class
    Foam::error

Description
    Class to handle errors and exceptions in a simple, consistent stream-based
    manner.

    The error class is globally instantiated with a title string. Errors,
    messages and other data are piped to the messageStream class in the
    standard manner.  Manipulators are supplied for exit and abort that may
    terminate the program or throw an exception depending on whether the
    exception handling has been switched on (off by default).

Usage
    \code
        error << "message1" << "message2" << FoamDataType << exit(errNo);
        error << "message1" << "message2" << FoamDataType << abort();
    \endcode

SourceFiles
    error.C
    IOerror.C

SeeAlso
    Foam::errorManip
    Foam::errorManipArg

\*---------------------------------------------------------------------------*/

#ifndef error_H
#define error_H

#include "messageStream.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class error Declaration
\*---------------------------------------------------------------------------*/

class error
:
    public std::exception,
    public messageStream
{
    // Private Member Functions

        //- Common code for exit or abort
        void exitOrAbort(const int errNo, const bool isAbort);


protected:

    // Protected Data

        string functionName_;
        string sourceFileName_;
        label sourceFileLineNumber_;
        bool throwing_;
        std::unique_ptr<OStringStream> messageStreamPtr_;


public:

    // Constructors

        //- Construct from title string
        explicit error(const string& title);

        //- Construct from dictionary
        explicit error(const dictionary& errDict);

        //- Copy construct
        error(const error& err);


    //- Destructor
    virtual ~error() noexcept;


    // Static Functions

        //- Emit warning on stderr about something being old.
        //  \param what description for the warning
        //  \param version is the old version (YYMM) for determining the
        //      age in months compared to the current OpenFOAM version
        //      as conveyed by the \c foamVersion::api value.
        static void warnAboutAge(const char* what, const int version);


    // Member Functions

        //- The accumulated error message
        string message() const;

        //- Clear any messages
        void clear() const;

        const string& functionName() const
        {
            return functionName_;
        }

        const string& sourceFileName() const
        {
            return sourceFileName_;
        }

        label sourceFileLineNumber() const
        {
            return sourceFileLineNumber_;
        }

        //- Return the current exception throwing (on or off)
        bool throwing() const
        {
            return throwing_;
        }

        //- Activate/deactivate exception throwing
        //  \return the previous throwing state
        inline bool throwExceptions(bool doThrow)
        {
            const bool prev = throwing_;
            throwing_ = doThrow;
            return prev;
        }

        //- Activate exception throwing
        //  \return the previous throwing state
        inline bool throwExceptions()
        {
            return throwExceptions(true);
        }

        //- Deactivate exception throwing
        //  \return the previous throwing state
        inline bool dontThrowExceptions()
        {
            return throwExceptions(false);
        }

        //- Define basic print message
        //  \return OSstream for further info.
        OSstream& operator()
        (
            const string& functionName
        );

        //- Define basic print message
        //  \return OSstream for further info.
        OSstream& operator()
        (
            const char* functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber = 0
        );

        //- Define basic print message
        //  \return OSstream for further info.
        OSstream& operator()
        (
            const string& functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber = 0
        );

        //- Explicit convert to OSstream for << operations
        OSstream& operator()()
        {
            return operator OSstream&();
        }

        //- Convert to OSstream
        operator OSstream&();

        //- Create and return a dictionary representation of the error
        operator dictionary() const;


        //- Helper function to print a stack,
        //- used when OpenFOAM IO not yet initialised.
        static void safePrintStack(std::ostream& os);

        //- Helper function to print a stack
        static void printStack(Ostream& os);

        //- True if FOAM_ABORT is on.
        static bool useAbort();

        //- Exit : can be called for any error to exit program.
        //  Redirects to abort() when FOAM_ABORT is on.
        void exit(const int errNo = 1);

        //- Abort : used to stop code for fatal errors.
        //  Prints stack before exiting.
        void abort();

        //- Print error message
        void write(Ostream& os, const bool includeTitle = true) const;
};


/*---------------------------------------------------------------------------*\
                           Class IOerror Declaration
\*---------------------------------------------------------------------------*/

//- Report an I/O error
class IOerror
:
    public error
{
    // Private Data

        string ioFileName_;
        label ioStartLineNumber_;
        label ioEndLineNumber_;


    // Private Member Functions

        //- Common code for exit or abort
        void exitOrAbort(const int errNo, const bool isAbort);


public:

    // Constructors

        //- Construct from title string
        explicit IOerror(const string& title);

        //- Construct from dictionary
        explicit IOerror(const dictionary& errDict);


    //- Destructor
    virtual ~IOerror() noexcept;


    // Member Functions

        const string& ioFileName() const
        {
            return ioFileName_;
        }

        label ioStartLineNumber() const
        {
            return ioStartLineNumber_;
        }

        label ioEndLineNumber() const
        {
            return ioEndLineNumber_;
        }

        //- Convert to OSstream
        //  Prints basic message and returns OSstream for further info.
        OSstream& operator()
        (
            const char* functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber,
            const string& ioFileName,
            const label ioStartLineNumber = -1,
            const label ioEndLineNumber = -1
        );

        //- Convert to OSstream
        //  Prints basic message and returns OSstream for further info.
        OSstream& operator()
        (
            const char* functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber,
            const IOstream& ioStream
        );

        //- Convert to OSstream
        //  Prints basic message and returns OSstream for further info.
        OSstream& operator()
        (
            const char* functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber,
            const dictionary& dict
        );

        //- Print basic message and exit.
        //  Uses cerr if streams not yet constructed (at startup).
        //  Use in startup parsing instead of FatalError.
        static void SafeFatalIOError
        (
            const char* functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber,
            const IOstream& ioStream,
            const string& msg
        );

        //- Create and return a dictionary representation of the error
        operator dictionary() const;


        //- Exit : can be called for any error to exit program
        void exit(const int errNo = 1);

        //- Abort : used to stop code for fatal errors
        void abort();

        //- Print error message
        void write(Ostream& os, const bool includeTitle = true) const;
};


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Ostream operator
Ostream& operator<<(Ostream& os, const error& err);

//- Ostream operator
Ostream& operator<<(Ostream& os, const IOerror& err);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global error declarations: defined in error.C

//- Error stream (uses stdout - output on all processes),
//- with additional 'FOAM FATAL ERROR' header text and stack trace.
extern error   FatalError;

//- Error stream (uses stdout - output on all processes),
//- with additional 'FOAM FATAL IO ERROR' header text and stack trace.
extern IOerror FatalIOError;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Convenience macros to add the file name and line number to the function name

//- Report an error message using Foam::FatalError
//  for functionName in file __FILE__ at line __LINE__
#define FatalErrorIn(functionName)                                             \
    ::Foam::FatalError((functionName), __FILE__, __LINE__)

//- Report an error message using Foam::FatalError
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
#define FatalErrorInFunction FatalErrorIn(FUNCTION_NAME)

//- Report an error message using Foam::FatalError
//  about unknown lookup type in table
#define FatalErrorInLookup(lookupTag, lookupName, lookupTable)                 \
    ::Foam::FatalError(FUNCTION_NAME, __FILE__, __LINE__)                      \
    << "Unknown " << (lookupTag) << " type " << (lookupName)                   \
    << "\n\nValid " << (lookupTag) << " types :\n"                             \
    << ((lookupTable).sortedToc()) << '\n'


//- Report an error message using Foam::FatalIOError
//  for functionName in file __FILE__ at line __LINE__
//  for a particular IOstream
#define FatalIOErrorIn(functionName, ios)                                      \
    ::Foam::FatalIOError((functionName), __FILE__, __LINE__, (ios))

//- Report an error message using Foam::FatalIOError
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
//  for a particular IOstream
#define FatalIOErrorInFunction(ios) FatalIOErrorIn(FUNCTION_NAME, ios)


//- Report an error message using Foam::FatalIOError
//  about unknown lookup type in table
#define FatalIOErrorInLookup(ios, lookupTag, lookupName, lookupTable)          \
    ::Foam::FatalIOError(FUNCTION_NAME, __FILE__, __LINE__, (ios))             \
    << "Unknown " << (lookupTag) << " type " << (lookupName)                   \
    << "\n\nValid " << (lookupTag) << " types :\n"                             \
    << ((lookupTable).sortedToc()) << '\n'


//- Report an error message using Foam::FatalIOError
//  (or cerr if FatalIOError not yet constructed)
//  for functionName in file __FILE__ at line __LINE__
//  for a particular IOstream
#define SafeFatalIOErrorIn(functionName, ios, msg)                             \
    ::Foam::IOerror::SafeFatalIOError                                          \
    ((functionName), __FILE__, __LINE__, (ios), (msg))

//- Report an error message using Foam::FatalIOError
//  (or cerr if FatalIOError not yet constructed)
//  for functionName in file __FILE__ at line __LINE__
//  for a particular IOstream
#define SafeFatalIOErrorInFunction(ios, msg)                                   \
    SafeFatalIOErrorIn(FUNCTION_NAME, ios, msg)


//- Issue a FatalErrorIn for a function not currently implemented.
//  The functionName is printed and then abort is called.
//
//  This macro can be particularly useful when methods must be defined to
//  complete the interface of a derived class even if they should never be
//  called for this derived class.
#define notImplemented(functionName)                                           \
    FatalErrorIn(functionName)                                                 \
        << "Not implemented" << ::Foam::abort(FatalError);

//- Issue a FatalErrorIn for a function not currently implemented.
//  The FUNCTION_NAME is printed and then abort is called.
//
//  This macro can be particularly useful when methods must be defined to
//  complete the interface of a derived class even if they should never be
//  called for this derived class.
#define NotImplemented notImplemented(FUNCTION_NAME)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "errorManip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
