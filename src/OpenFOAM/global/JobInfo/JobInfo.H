/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::JobInfo

Description
    Helper class for recording information about run/finished jobs.

    Writes the following files:
      - $FOAM_JOB_DIR/runningJobs
      - $FOAM_JOB_DIR/finishedJobs

    If FOAM_JOB_DIR is unset, defaults to ~/.OpenFOAM/jobControl

SourceFiles
    JobInfo.C

\*---------------------------------------------------------------------------*/

#ifndef JobInfo_H
#define JobInfo_H

#include "dictionary.H"
#include "fileName.H"
#include "cpuTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class JobInfo Declaration
\*---------------------------------------------------------------------------*/

class JobInfo
:
    public dictionary
{
    // Private data

        //- The name of the job file
        string   jobFileName_;
        fileName runningDir_;
        fileName finishedDir_;
        cpuTime  cpuTime_;

    // Private Member Functions

        bool write(Ostream& os) const;
        void end(const word& terminationType);


public:

    //- Global value for constructed job info
    static bool constructed;

    //- Global value for writeJobInfo enabled
    static bool writeJobInfo;


    // Constructors

        //- Construct null
        JobInfo();


    //- Destructor
    //  Update job info and relocate the file from running to finished.
    ~JobInfo();


    // Member Functions

        //- Write the job info to its file in the runningJobs directory
        void write() const;

        //- End with "termination=normal"
        void end();

        //- End with "termination=exit"
        void exit();

        //- End with "termination=abort"
        void abort();

        //- Update job info and relocate the file from running to finished.
        void signalEnd() const;
};


extern JobInfo jobInfo;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
