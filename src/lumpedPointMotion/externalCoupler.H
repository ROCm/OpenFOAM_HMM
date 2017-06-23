/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify i
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
    Foam::externalCoupler

Description
    Encapsulates the logic for coordinating between OpenFOAM and an
    external application.

    This class provides a simple interface for explicit coupling with an
    external application using plain text files located in the user-specified
    communications directory.
    These files are to be read/written on the master processor only.

    The explicit coupling follows a master/slave model in which OpenFOAM
    is the 'master' and the external application is the 'slave'.
    The readiness to exchange information in either direction is handled
    by a lock file (ie, OpenFOAM.lock).
    If the lock file is present, the slave (external application) should wait
    for the master (OpenFOAM) to complete.

    When the master is finished its tasks and has prepared data for
    the slave, the lock file is removed, instructing the external
    source to take control of the program execution.

    When the slave has completed its tasks, it will reinstate the lock file.

    Example of the communication specification:
    \verbatim
    communication
    {
        commsDir        "${FOAM_CASE}/comms";
        waitInterval    1;
        timeOut         100;
        initByExternal  no;
    }
    \endverbatim

SourceFiles
    externalCoupler.C

\*---------------------------------------------------------------------------*/

#ifndef externalCoupler_H
#define externalCoupler_H

#include "fileName.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class externalCoupler Declaration
\*---------------------------------------------------------------------------*/

class externalCoupler
{
    //- The run state (ie, who is currently in charge)
    enum runState
    {
        NONE,       //!< Not initialized
        MASTER,     //!< The master (OpenFOAM) is in charge
        SLAVE,      //!< The slave (external program) is in charge
        DONE        //!< Finished
    };


    // Private data

        //- The current run (and initialization) state
        mutable runState runState_;

        //- Local path to communications directory
        fileName commsDir_;

        //- Interval time between checking for return data [s]
        unsigned waitInterval_;

        //- Timeout [s] while waiting for the external application
        unsigned timeOut_;

        //- Flag to indicate values are initialized by external application
        bool slaveFirst_;

        //- Local logging/verbosity flag
        bool log;


    // Private Member Functions

        //- Return the file path to the base communications directory
        const fileName& baseDir() const;

        //- Return the file path to the lock file
        fileName lockFile() const;

        //- Disallow default bitwise copy construc
        externalCoupler(const externalCoupler&) = delete;

        //- Disallow default bitwise assignmen
        void operator=(const externalCoupler&) = delete;

public:

    //- Runtime type information
    TypeName("externalCoupler");


    // Static data members

        //- Name of the lock file
        static word lockName;


    // Constructors

        //- Construct null using standard defaults
        externalCoupler();

        //- Construct from dictionary
        externalCoupler(const dictionary& dict);


    //- Destructor
    virtual ~externalCoupler();


    // Member Functions

      // Access

        //- True if state has been initialized
        bool initialized() const;

        //- External application provides initial values
        bool slaveFirst() const;

        //- Create lock file to indicate that OpenFOAM is in charge
        //  Optionally wait for master as well.
        void useMaster(const bool wait=false) const;

        //- Wait for indication that OpenFOAM has supplied output.
        //  This is when the lock file disappears, or it exists but with
        //  "status=done" content.
        //  \return False if lock file contains "status=done"
        bool waitForMaster() const;

        //- Remove lock file to indicate that the external program is in charge
        //  Optionally wait for slave as well.
        void useSlave(const bool wait=false) const;

        //- Wait for indication that the external program has supplied input.
        //  This is when the lock file appears.
        //  \return False if lock file contains "status=done"
        bool waitForSlave() const;

        //- Return the file path in the communications directory
        fileName resolveFile(const word& file) const;

        //- Generate status=done in lock (only when run-state = master)
        void shutdown() const;

        //- Remove files written by OpenFOAM
        void removeDirectory() const;


      // Edit

        //- Read communication settings from dictionary
        bool readDict(const dictionary& dict);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
