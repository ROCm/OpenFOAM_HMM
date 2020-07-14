/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::Time

Description
    Class to control time during OpenFOAM simulations that is also the
    top-level objectRegistry.

SourceFiles
    Time.C
    TimeIO.C
    findInstance.C

\*---------------------------------------------------------------------------*/

#ifndef Time_H
#define Time_H

#include "TimePaths.H"
#include "objectRegistry.H"
#include "unwatchedIOdictionary.H"
#include "FIFOStack.H"
#include "clock.H"
#include "cpuTime.H"
#include "TimeState.H"
#include "Switch.H"
#include "instantList.H"
#include "Enum.H"
#include "typeInfo.H"
#include "dlLibraryTable.H"
#include "functionObjectList.H"
#include "sigWriteNow.H"
#include "sigStopAtWriteNow.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class argList;
class profilingTrigger;
class OSstream;

/*---------------------------------------------------------------------------*\
                             Class Time Declaration
\*---------------------------------------------------------------------------*/

class Time
:
    public clock,
    public cpuTime,
    public TimePaths,
    public objectRegistry,
    public TimeState
{
public:

        //- Write control options
        enum writeControls
        {
            wcNone,                 //!< "none"
            wcTimeStep,             //!< "timeStep"
            wcRunTime,              //!< "runTime"
            wcAdjustableRunTime,    //!< "adjustable" / "adjustableRunTime"
            wcClockTime,            //!< "clockTime"
            wcCpuTime,              //!< "cpuTime"
            wcUnknown               //!< Dummy no-op
        };

        //- Stop-run control options, which are primarily used when
        //- altering the stopAt condition.
        enum stopAtControls
        {
            saEndTime,    //!< Stop when Time reaches prescribed endTime
            saNoWriteNow, //!< Adjust endTime to stop immediately w/o writing
            saWriteNow,   //!< adjust endTime to stop immediately w/ writing
            saNextWrite,  //!< stop at the next data write interval
            saUnknown     //!< Dummy no-op. Do not change current value.
        };

        //- Supported time directory name formats
        enum fmtflags
        {
            general    = 0,                   //!< default float notation
            fixed      = ios_base::fixed,     //!< fixed-point notation
            scientific = ios_base::scientific //!< scientific notation
        };


        //- Names for writeControls
        static const Enum<writeControls> writeControlNames;

        //- Names for stopAtControls
        static const Enum<stopAtControls> stopAtControlNames;


        //- Style for "ExecutionTime = " output
        //  0 = seconds (with trailing 's')
        //  1 = day-hh:mm:ss
        //
        // \note this is public so registered info switches can modify it.
        static int printExecutionFormat_;


private:

    // Private Data

        //- Profiling trigger for time-loop (for run, loop)
        mutable std::unique_ptr<profilingTrigger> loopProfiling_;

        //- Any loaded dynamic libraries
        //  Construct before reading controlDict
        mutable dlLibraryTable libs_;

        //- The controlDict
        unwatchedIOdictionary controlDict_;


protected:

    // Protected Data

        label startTimeIndex_;

        scalar startTime_;

        mutable scalar endTime_;

        mutable stopAtControls stopAt_;

        writeControls writeControl_;

        scalar writeInterval_;

        label purgeWrite_;

        mutable FIFOStack<word> previousWriteTimes_;

        //- The total number of sub-cycles, the current sub-cycle index,
        //- or 0 if time is not being sub-cycled
        label subCycling_;

        // One-shot writing
        bool writeOnce_;

        //- If time is being sub-cycled this is the previous TimeState
        autoPtr<TimeState> prevTimeState_;

        //- Signal handler for one-shot writing upon signal
        sigWriteNow sigWriteNow_;

        //- Signal handler for write and clean exit upon signal
        sigStopAtWriteNow sigStopAtWriteNow_;


        //- Time directory name format
        static fmtflags format_;

        //- Time directory name precision
        static int precision_;

        //- Maximum time directory name precision
        static const int maxPrecision_;


        //- Adjust the time step so that writing occurs at the specified time
        void adjustDeltaT();

        //- Set the controls from the current controlDict
        void setControls();

        //- Set file monitoring, profiling, etc
        //  Optionally force profiling without inspecting the controlDict
        void setMonitoring(const bool forceProfiling=false);

        //- Read the control dictionary and set the write controls etc.
        virtual void readDict();


private:

        //- The write stream option (format, compression, version)
        IOstreamOption writeStreamOption_;

        //- Default graph format
        word graphFormat_;

        //- Is runtime modification of dictionaries allowed?
        Switch runTimeModifiable_;

        //- Function objects executed at start and on ++, +=
        mutable functionObjectList functionObjects_;


public:

    // Static Data Members

        //- The default control dictionary name (normally "controlDict")
        static word controlDictName;


    //- Runtime type information
    TypeName("time");


    // Constructors

        //- Construct from dictionary name to read and argument list.
        //  Uses "system" and "constant" for its directories.
        inline Time
        (
            const word& ctrlDictName,
            const argList& args,
            const bool enableFunctionObjects = true,
            const bool enableLibs = true
        );

        //- Construct given name of dictionary to read and root/case paths.
        //  Uses "system" and "constant" for its directories.
        inline Time
        (
            const word& ctrlDictName,
            const fileName& rootPath,
            const fileName& caseName,
            const bool enableFunctionObjects = true,
            const bool enableLibs = true
        );

        //- Construct given dictionary and root/case paths.
        //  Uses "system" and "constant" for its directories.
        inline Time
        (
            const dictionary& dict,
            const fileName& rootPath,
            const fileName& caseName,
            const bool enableFunctionObjects = true,
            const bool enableLibs = true
        );

        //- Construct at endTime and root/case paths.
        //  Uses "system" and "constant" for its directories.
        inline Time
        (
            const fileName& rootPath,
            const fileName& caseName,
            const bool enableFunctionObjects = true,
            const bool enableLibs = true
        );

        //- Construct from dictionary name to read and argument list.
        //  Use specified system/constant directory names.
        Time
        (
            const word& ctrlDictName,
            const argList& args,
            const word& systemName,     //!< "system" directory name
            const word& constantName,   //!< "constant" directory name
            const bool enableFunctionObjects = true,
            const bool enableLibs = true
        );

        //- Construct given name of dictionary to read and root/case paths.
        //  Use specified system/constant directory names.
        Time
        (
            const word& ctrlDictName,
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName,     //!< "system" directory name
            const word& constantName,   //!< "constant" directory name
            const bool enableFunctionObjects = true,
            const bool enableLibs = true
        );

        //- Construct given dictionary and root/case paths.
        //  Use specified system/constant directory names.
        Time
        (
            const dictionary& dict,
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName,     //!< "system" directory name
            const word& constantName,   //!< "constant" directory name
            const bool enableFunctionObjects = true,
            const bool enableLibs = true
        );

        //- Construct given endTime and root/case paths.
        //  Use specified system/constant directory names.
        Time
        (
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName,     //!< "system" directory name
            const word& constantName,   //!< "constant" directory name
            const bool enableFunctionObjects = true,
            const bool enableLibs = true
        );


    // Selectors

        //- Construct (dummy) Time - no functionObjects or libraries
        static autoPtr<Time> New();

        //- Construct (dummy) Time - no functionObjects or libraries
        static autoPtr<Time> New(const fileName& caseDir);

        //- Construct (dummy) Time - no functionObjects or libraries
        static autoPtr<Time> New(const argList& args);


    //- Destructor
    virtual ~Time();


    // Member Functions

        // Database functions

            //- Return name from objectRegistry and not TimePaths
            using objectRegistry::name;

            //- Return root path
            using TimePaths::rootPath;

            //- Return case name
            using TimePaths::caseName;

            //- Return path
            fileName path() const
            {
                return rootPath()/caseName();
            }

            //- Return read access to the controlDict dictionary
            const dictionary& controlDict() const
            {
                return controlDict_;
            }

            virtual const fileName& dbDir() const
            {
                return fileName::null;
            }

            //- Return current time path
            fileName timePath() const
            {
                return path()/timeName();
            }

            //- The write stream option (format, compression, version)
            IOstreamOption writeStreamOption() const
            {
                return writeStreamOption_;
            }

            //- The write stream format
            IOstream::streamFormat writeFormat() const
            {
                return writeStreamOption_.format();
            }

            //- The write stream compression
            IOstream::compressionType writeCompression() const
            {
                return writeStreamOption_.compression();
            }

            //- The write stream version
            IOstream::versionNumber writeVersion() const
            {
                return writeStreamOption_.version();
            }

            //- Default graph format
            const word& graphFormat() const
            {
                return graphFormat_;
            }

            //- Supports re-reading
            const Switch& runTimeModifiable() const
            {
                return runTimeModifiable_;
            }

            //- Read control dictionary, update controls and time
            virtual bool read();

            //- Read the objects that have been modified
            void readModifiedObjects();

            //- Return the location of "dir" containing the file "name".
            //- (eg, used in reading mesh data)
            //  If name is null, search for the directory "dir" only.
            //  Does not search beyond stopInstance (if set) or constant.
            word findInstance
            (
                const fileName& dir,
                const word& name = word::null,
                const IOobject::readOption rOpt = IOobject::MUST_READ,
                const word& stopInstance = word::null
            ) const;

            //- Search the case for the time directory path
            //- corresponding to the given instance
            word findInstancePath
            (
                const fileName& directory,
                const instant& t
            ) const;

            //- Search the case for the time directory path
            //- corresponding to the given instance
            word findInstancePath(const instant& t) const;

            //- Write time dictionary to the \<time\>/uniform directory
            virtual bool writeTimeDict() const;

            //- Write using stream options
            virtual bool writeObject
            (
                IOstreamOption streamOpt,
                const bool valid
            ) const;

            //- Write the objects immediately (not at end of iteration)
            //- and continue the run
            bool writeNow();

            //- Write the objects now (not at end of iteration) and end the run
            bool writeAndEnd();

            //- Write the objects once (one shot) and continue the run
            void writeOnce();

            //- Print the elapsed ExecutionTime (cpu-time), ClockTime
            Ostream& printExecutionTime(OSstream& os) const;


        // Access

            //- Return time name of given scalar time
            //- formatted with the given precision
            static word timeName
            (
                const scalar t,
                const int precision = precision_
            );

            //- Return current time name
            virtual word timeName() const;

            //- Return start time index
            virtual label startTimeIndex() const;

            //- Return start time
            virtual dimensionedScalar startTime() const;

            //- Return end time
            virtual dimensionedScalar endTime() const;

            //- Return the stop control information
            virtual stopAtControls stopAt() const;

            //- Return true if adjustTimeStep is true
            virtual bool isAdjustTimeStep() const;

            //- Return the list of function objects
            const functionObjectList& functionObjects() const
            {
                return functionObjects_;
            }

            //- Mutable access to the loaded dynamic libraries
            dlLibraryTable& libs() const
            {
                return libs_;
            }

            //- Zero (tests as false) if time is not being sub-cycled,
            //- otherwise the current sub-cycle index or the total number of
            //- sub-cycles.
            //  The interpretation of non-zero values is dependent on the
            //  routine.
            label subCycling() const
            {
                return subCycling_;
            }

            //- Return previous TimeState if time is being sub-cycled
            const TimeState& prevTimeState() const
            {
                return *prevTimeState_;
            }


        // Check

            //- Return true if run should continue,
            //  also invokes the functionObjectList::end() method
            //  when the time goes out of range
            //  \note
            //  For correct behaviour, the following style of time-loop
            //  is recommended:
            //  \code
            //      while (runTime.run())
            //      {
            //          ++runTime;
            //          solve;
            //          runTime.write();
            //      }
            //  \endcode
            virtual bool run() const;

            //- Return true if run should continue and if so increment time
            //  also invokes the functionObjectList::end() method
            //  when the time goes out of range
            //  \note
            //  For correct behaviour, the following style of time-loop
            //  is recommended:
            //  \code
            //      while (runTime.loop())
            //      {
            //          solve;
            //          runTime.write();
            //      }
            //  \endcode
            virtual bool loop();

            //- Return true if end of run,
            //  does not invoke any functionObject methods
            //  \note
            //      The rounding heuristics near endTime mean that
            //      \code run() \endcode and \code !end() \endcode may
            //      not yield the same result
            virtual bool end() const;


        // Edit

            //- Adjust the current stopAtControl.
            //  \param stopCtrl the new stop control, whereby
            //      stopAtControls::saUnknown is treated as a no-op.
            //  \note this value only persists until the next time the
            //      dictionary is read.
            //  \return true if the stopAt() value was changed.
            virtual bool stopAt(const stopAtControls stopCtrl) const;

            //- Reset the time and time-index to those of the given time
            virtual void setTime(const Time& t);

            //- Reset the time and time-index
            virtual void setTime(const instant& inst, const label newIndex);

            //- Reset the time and time-index
            virtual void setTime
            (
                const dimensionedScalar& newTime,
                const label newIndex
            );

            //- Reset the time and time-index
            virtual void setTime(const scalar newTime, const label newIndex);

            //- Reset end time
            virtual void setEndTime(const dimensionedScalar& endTime);

            //- Reset end time
            virtual void setEndTime(const scalar endTime);

            //- Reset time step, normally also calling adjustDeltaT()
            virtual void setDeltaT
            (
                const dimensionedScalar& deltaT,
                const bool adjust = true
            );

            //- Reset time step, normally also calling adjustDeltaT()
            virtual void setDeltaT
            (
                const scalar deltaT,
                const bool adjust = true
            );

            //- Set time to sub-cycle for the given number of steps
            virtual TimeState subCycle(const label nSubCycles);

            //- Adjust the reported sub-cycle index.
            //  \param index is the sub-cycle index.
            //      This index is ignored sub-cycling was
            //      not already registered, or if the index is zero or
            //      negative.
            virtual void subCycleIndex(const label index);

            //- Reset time after sub-cycling back to previous TimeState
            virtual void endSubCycle();

            //- Return non-const access to the list of function objects
            functionObjectList& functionObjects()
            {
                return functionObjects_;
            }


    // Member operators

        //- Set deltaT to that specified and increment time via operator++()
        virtual Time& operator+=(const dimensionedScalar& deltaT);

        //- Set deltaT to that specified and increment time via operator++()
        virtual Time& operator+=(const scalar deltaT);

        //- Prefix increment,
        //  also invokes the functionObjectList::start() or
        //  functionObjectList::execute() method, depending on the time-index
        virtual Time& operator++();

        //- Postfix increment, this is identical to the prefix increment
        virtual Time& operator++(int);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TimeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
