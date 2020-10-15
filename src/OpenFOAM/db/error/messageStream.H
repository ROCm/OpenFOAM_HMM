/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::messageStream

Description
    Class to handle messaging in a simple, consistent stream-based
    manner.

    The messageStream class is globally instantiated with a title string and
    a severity (which controls the program termination) and a number of
    errors before termination.  Errors, messages and other data are sent to
    the messageStream class in the standard manner.

Usage
    \code
        messageStream
            << "message1" << "message2" << FoamDataType << endl;
    \endcode

SourceFiles
    messageStream.C

\*---------------------------------------------------------------------------*/

#ifndef messageStream_H
#define messageStream_H

#include "label.H"
#include "string.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class IOstream;
class Ostream;
class OSstream;
class OStringStream;
class dictionary;

/*---------------------------------------------------------------------------*\
                           Class messageStream Declaration
\*---------------------------------------------------------------------------*/

class messageStream
{
public:

    //- Message type, or error severity flags
    enum errorSeverity
    {
        INFO = 1,   //!< General information output
        WARNING,    //!< Warning of possible problem.
        SERIOUS,    //!< A serious problem - eg, data corruption.
        FATAL,      //!< A fatal error.
        INFO_STDERR = INFO | 0x10,  //!< Information, but on stderr
    };


protected:

    // Protected Data

        string title_;
        errorSeverity severity_;
        int maxErrors_;
        int errorCount_;


public:

    // Static Data

        //- The output level (verbosity) of messages
        //
        //  - level == 0 : suppress all output
        //  - level == 1 : normal output
        //  - level >= 2 : report source file name and line number if available
        //
        // \note The default level is normally 2.
        static int level;

        //- The output redirection of messages
        //
        //  - redirect == 2 : use stderr instead of stdout
        static int redirect;


    // Constructors

        //- Construct from components
        messageStream
        (
            const string& title,
            const errorSeverity severity,
            const int maxErrors = 0
        );

        //- Construct as Fatal from dictionary, extracting the 'title'.
        explicit messageStream(const dictionary& dict);


    // Member functions

        //- The title of this error type
        const string& title() const
        {
            return title_;
        }

        //- The maximum number of errors before program termination
        int maxErrors() const
        {
            return maxErrors_;
        }

        //- Non-const access to the maximum number of errors before
        //- program termination to enable user to reset it
        int& maxErrors()
        {
            return maxErrors_;
        }

        //- Convert to OSstream
        //  Prints to Pout for the master stream
        OSstream& masterStream(const label communicator);

        //- Print basic message
        //  \return OSstream for further info.
        OSstream& operator()
        (
            const string& functionName
        );

        //- Print basic message
        //  \return OSstream for further info.
        OSstream& operator()
        (
            const char* functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber = 0
        );

        //- Print basic message
        //  \return OSstream for further info.
        OSstream& operator()
        (
            const string& functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber = 0
        );

        //- Print basic message
        //  \return OSstream for further info.
        OSstream& operator()
        (
            const char* functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber,
            const string& ioFileName,
            const label ioStartLineNumber = -1,
            const label ioEndLineNumber = -1
        );

        //- Print basic message
        //  \return OSstream for further info.
        OSstream& operator()
        (
            const char* functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber,
            const IOstream&
        );

        //- Print basic message
        //  \return OSstream for further info.
        OSstream& operator()
        (
            const char* functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber,
            const dictionary&
        );

        //- Convert to OSstream for << operations
        operator OSstream&();

        //- Explicitly convert to OSstream for << operations
        OSstream& operator()()
        {
            return operator OSstream&();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global error declarations: defined in messageStream.C

//- Global for selective suppression of Info output.
//  This is normally accessed implicitly via the DetailInfo macro and is often
//  associated applications with suppressed banners. For example,
//
//  \code
//      DetailInfo << "Hello, I'm running from program xyz" << nl;
//      Info<< "Found ... invalid items" << nl;
//  \endcode
//
//  The values are normally 0 or a positive value.
//  \note This flag is initialized to 1 by default.
extern int infoDetailLevel;

//- Information stream (uses stdout - output is on the master only)
extern messageStream Info;

//- Information stream (uses stderr - output is on the master only)
extern messageStream InfoErr;

//- Warning stream (uses stdout - output is on the master only),
//- with additional 'FOAM Warning' header text.
extern messageStream Warning;

//- Error stream (uses stdout - output on all processes),
//- with additional 'FOAM Serious Error' header text.
extern messageStream SeriousError;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "OSstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Convenience macros to add the file name and line number to the function name

// Compiler provided function name string:
//     for gcc-compatible compilers use __PRETTY_FUNCTION__
//     otherwise use the standard __func__
#ifdef __GNUC__
    #define FUNCTION_NAME __PRETTY_FUNCTION__
#else
    #define FUNCTION_NAME __func__
#endif


//- Report an error message using Foam::SeriousError
//  for functionName in file __FILE__ at line __LINE__
#define SeriousErrorIn(functionName)                                           \
    ::Foam::SeriousError((functionName), __FILE__, __LINE__)

//- Report an error message using Foam::SeriousError
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
#define SeriousErrorInFunction SeriousErrorIn(FUNCTION_NAME)


//- Report an IO error message using Foam::SeriousError
//  for functionName in file __FILE__ at line __LINE__
//  for a particular IOstream
#define SeriousIOErrorIn(functionName, ios)                                    \
    ::Foam::SeriousError((functionName), __FILE__, __LINE__, ios)

//- Report an IO error message using Foam::SeriousError
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
//  for a particular IOstream
#define SeriousIOErrorInFunction(ios) SeriousIOErrorIn(FUNCTION_NAME, ios)


//- Report a warning using Foam::Warning
//  for functionName in file __FILE__ at line __LINE__
#define WarningIn(functionName)                                                \
    ::Foam::Warning((functionName), __FILE__, __LINE__)

//- Report a warning using Foam::Warning
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
#define WarningInFunction WarningIn(FUNCTION_NAME)


//- Report an IO warning using Foam::Warning
//  for functionName in file __FILE__ at line __LINE__
//  for a particular IOstream
#define IOWarningIn(functionName, ios)                                         \
    ::Foam::Warning((functionName), __FILE__, __LINE__, (ios))

//- Report an IO warning using Foam::Warning
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
//  for a particular IOstream
#define IOWarningInFunction(ios) IOWarningIn(FUNCTION_NAME, ios)


//- Report an information message using Foam::Info
//  for functionName in file __FILE__ at line __LINE__
#define InfoIn(functionName)                                                   \
    ::Foam::Info((functionName), __FILE__, __LINE__)

//- Report an information message using Foam::Info
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
#define InfoInFunction InfoIn(FUNCTION_NAME)

//- Report an information message using Foam::Pout
//  for functionName in file __FILE__ at line __LINE__
#define PoutIn(functionName)                                                   \
    ::Foam::Pout((functionName), __FILE__, __LINE__)

//- Report an information message using Foam::Pout
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
#define PoutInFunction PoutIn(FUNCTION_NAME)

//- Write to Foam::Info if the Foam::infoDetailLevel is +ve non-zero (default)
#define DetailInfo                                                             \
    if (::Foam::infoDetailLevel > 0) Info

//- Report write to Foam::Info if the local log switch is true
#define Log                                                                    \
    if (log) Info


//- Report an IO information message using Foam::Info
//  for functionName in file __FILE__ at line __LINE__
//  for a particular IOstream
#define IOInfoIn(functionName, ios)                                            \
    ::Foam::Info((functionName), __FILE__, __LINE__, (ios))

//- Report an IO information message using Foam::Info
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
//  for a particular IOstream
#define IOInfoInFunction(ios) IOInfoIn(FUNCTION_NAME, ios)


//- Report an information message using Foam::Info
//  if the local debug switch is true
#define DebugInfo                                                              \
    if (debug) Info

//- Report an information message using Foam::Info
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
//  if the local debug switch is true
#define DebugInFunction                                                        \
    if (debug) InfoInFunction

//- Report an information message using Foam::Pout
//  if the local debug switch is true
#define DebugPout                                                              \
    if (debug) Pout

//- Report an information message using Foam::Pout
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
//  if the local debug switch is true
#define DebugPoutInFunction                                                    \
    if (debug) PoutInFunction

//- Report a variable name and value
//  using Foam::Pout in file __FILE__ at line __LINE__
#define DebugVar(var)                                                          \
{                                                                              \
    ::Foam::string oldPrefix(::Foam::Pout.prefix());                           \
    ::Foam::Pout<< "["<< __FILE__ << ":" << __LINE__ << "] ";                  \
    ::Foam::Pout.prefix() = oldPrefix + #var " ";                              \
    ::Foam::Pout<< var << ::Foam::endl;                                        \
    ::Foam::Pout.prefix() = oldPrefix;                                         \
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
