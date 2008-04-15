/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Class
    Foam::ProcessControl

Description

SourceFiles
    ProcessControl.C

\*---------------------------------------------------------------------------*/

#ifndef ProcessControl_H
#define ProcessControl_H

#include "label.H"
#include "stringList.H"
#include "word.H"
#include "fileName.H"

#include <sys/types.h>
#include <signal.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

/*---------------------------------------------------------------------------*\
                           Class ProcessControl Declaration
\*---------------------------------------------------------------------------*/

class ProcessControl
{
    // Private data

        Foam::fileName remoteShell_;
        Foam::fileName remoteCp_;
        Foam::label timeOut_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ProcessControl(const ProcessControl&);

        //- Disallow default bitwise assignment
        void operator=(const ProcessControl&);

        //- Read the entries from the given dictionary if they are present
        void readEntries(const Foam::dictionary& procDict);

        static void sigAlarm(const int);


public:

    // Public data types

        //- error value for time outs
        static const int TIMEDOUT = -2;

        class sigAlarmException
        {};


    // Constructors

        //- Construct from configuration file
        ProcessControl
        (
            const Foam::fileName& fxSystemConfigFileName
        );


    // Member Functions

        // Access

            const Foam::fileName& remoteShell() const
            {
                return remoteShell_;
            }

            const Foam::fileName& remoteCp() const
            {
                return remoteCp_;
            }

            Foam::label timeOut() const
            {
                return timeOut_;
            }

            Foam::stringList remoteShellArgs
            (
                const Foam::string& userName,
                const Foam::string& hostName,
                const Foam::stringList& arguments,
                const Foam::string& logName,
                const bool backGround
            ) const;

            Foam::stringList remoteCpArgs
            (
                const Foam::string& userName,
                const Foam::string& hostName,
                const Foam::fileName& src,
                const Foam::fileName& dest
            ) const;


        //- Concaternate arList into a single string
        static Foam::string commandString(const Foam::stringList& argList);

        //- Fork command
        static pid_t fork
        (
            const Foam::stringList& argList,
            const Foam::fileName& logFile
        );

        //- Waitpid command
        static pid_t waitpid(pid_t pid);

        //- Run command using system
        static int system
        (
            const Foam::stringList& argList,
            const int timeOut=0
        );

        //- Kill job with given pid
        static int kill(pid_t pid, int sig=SIGINT);

        //- Kill job on given machine with given pid
        int kill(const Foam::word& host, pid_t pid, int sig=SIGINT) const;

        //- Suspend job on given machine with given pid
        int suspend(const Foam::word& host, pid_t pid) const;

        //- Continue job on given machine with given pid
        int cont(const Foam::word& host, pid_t pid) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
