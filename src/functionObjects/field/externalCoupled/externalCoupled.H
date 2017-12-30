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
    Foam::functionObjects::externalCoupled

Group
    grpFieldFunctionObjects

Description
    This functionObject provides a simple interface for explicit coupling with
    an external application.

    The coupling is through plain text files where OpenFOAM boundary data
    is read/written as one line per face (data from all processors collated):
    \verbatim
        # Patch: <patch name>
        <fld1> <fld2> .. <fldn>             //face0
        <fld1> <fld2> .. <fldn>             //face1
        ..
        <fld1> <fld2> .. <fldn>             //faceN
    \endverbatim

    where the actual entries depend on the bc type:
    - mixed: value, snGrad, refValue, refGrad, valueFraction
    - externalCoupledMixed: output of writeDataMaster
    - other: value, snGrad

    These text files are located in a user specified communications directory
    which gets read/written on the master processor only. In the
    communications directory the structure will be
    \verbatim
        <regionsName>/<patchGroup>/<fieldName>.[in|out]
    \endverbatim

    (where regionsName is either the name of a single region or a composite
    of multiple region names)

    At start-up, the boundary creates a lock file, i.e..
    \verbatim
        OpenFOAM.lock
    \endverbatim

    ... to signal the external source to wait. During the functionObject
    execution the boundary values are written to files (one per region,
    per patch(group), per field), e.g.
    \verbatim
        <regionsName>/<patchGroup>/<fieldName>.out
    \endverbatim

    The lock file is then removed, instructing the external source to take
    control of the program execution. When ready, the external program
    should create the return values, e.g. to files
    \verbatim
        <regionsName>/<patchGroup>/<fieldName>.in
    \endverbatim

    ... and then reinstate the lock file. The functionObject will then
    read these values, apply them to the boundary conditions and pass
    program execution back to OpenFOAM.

Usage
    \verbatim
    externalCoupled
    {
        type            externalCoupled;
        ...
        log             yes;
        commsDir        "${FOAM_CASE}/comms";
        initByExternal  yes;
        stateEnd        remove;  // (remove | done)

        regions
        {
            "(region1|region0)"         // Name of region(s)
            {
                TPatchGroup             // Name of patch(group)
                {
                    readFields  (p);    // List of fields to read
                    writeFields (T);    // List of fields to write
                }
            }
        }
    }
    \endverbatim


    This reads/writes (on the master processor) the directory:
    \verbatim
        comms/region0_region1/TPatchGroup/
    \endverbatim

    with contents:
    \verbatim
        patchPoints     (collected points)
        patchFaces      (collected faces)
        p.in            (input file of p, written by external application)
        T.out           (output file of T, written by OpenFOAM)
    \endverbatim

    The patchPoints/patchFaces files denote the (collated) geometry
    which will be written if it does not exist yet or can be written as
    a preprocessing step using the createExternalCoupledPatchGeometry
    application.

    The entries comprise:
    \table
        Property     | Description             | Required | Default value
        type         | type name: externalCoupled | yes |
        commsDir     | communication directory | yes |
        waitInterval | wait interval in (s)    | no | 1
        timeOut      | timeout in (s)          | no | 100*waitInterval
        stateEnd     | Lockfile treatment on termination | no | done
        initByExternal | initialization values supplied by external app | yes
        calcFrequency | calculation frequency  | no | 1
        regions  | the regions to couple | yes
    \endtable


SourceFiles
    externalCoupled.C
    externalCoupledTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_externalCoupled_H
#define functionObjects_externalCoupled_H

#include "functionObject.H"
#include "externalFileCoupler.H"
#include "DynamicList.H"
#include "wordReList.H"
#include "scalarField.H"
#include "Enum.H"
#include "Switch.H"
#include "UPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;
class IFstream;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
              Class functionObjects::externalCoupled Declaration
\*---------------------------------------------------------------------------*/

class externalCoupled
:
    public functionObject,
    public externalFileCoupler
{
public:

    // Public data types

        //- Lockfile state on termination
        enum stateEnd
        {
            REMOVE,     //!< Remove lock file on end
            DONE,       //!< Lock file contains status=done on end
            IGNORE      //!< Internal use only (for handling cleanup).
        };

private:

        //- State end names (NB, only selectable values itemized)
        static const Enum<stateEnd> stateEndNames_;

    // Private data

        //- Reference to the time database
        const Time& time_;

        //- Calculation frequency
        label calcFrequency_;

        //- The last timeIndex when coupling was triggered
        label lastTrigger_;

        //- Names of (composite) regions
        DynamicList<word> regionGroupNames_;

        // Per (composite) region the names of the regions
        DynamicList<wordList> regionGroupRegions_;

        // Per (composite) region the indices of the group information
        HashTable<labelList> regionToGroups_;

        // Per group the names of the patches/patchGroups
        DynamicList<wordRe> groupNames_;

        // Per group the names of the fields to read
        DynamicList<wordList> groupReadFields_;

        // Per group the names of the fields to write
        DynamicList<wordList> groupWriteFields_;

        //- Initialised coupling
        bool initialisedCoupling_;


    // Private Member Functions

        //- Return the file path to the communications directory for the region
        static fileName groupDir
        (
            const fileName& commsDir,
            const word& regionsName,
            const wordRe& groupName
        );

        //- Read data for a single region, single field
        template<class Type>
        bool readData
        (
            const UPtrList<const fvMesh>& meshes,
            const wordRe& groupName,
            const word& fieldName
        );

        //- Write data for a single region, single field
        template<class Type>
        bool writeData
        (
            const UPtrList<const fvMesh>& meshes,
            const wordRe& groupName,
            const word& fieldName
        ) const;

        void initCoupling();

        //- Read (and distribute) scalar columns from stream. Every processor
        //  gets nRows (= patch size) of these. Note: could make its argument
        //  ISstream& but then would need additional logic to construct valid
        //  stream on all processors.
        void readColumns
        (
            const label nRows,
            const label nColumns,
            autoPtr<IFstream>& masterFilePtr,
            List<scalarField>& data
        ) const;

        //- Read (and distribute) lines from stream. Every processor
        //  gets nRows (= patch size) of these. Data kept as stream (instead
        //  of strings) for ease of interfacing to readData routines that take
        //  an Istream.
        void readLines
        (
            const label nRows,
            autoPtr<IFstream>& masterFilePtr,
            OStringStream& data
        ) const;

        //- Helper: append data from all processors onto master
        template<class Type>
        static tmp<Field<Type>> gatherAndCombine(const Field<Type>& fld);

        static void checkOrder(const wordList& regionNames);

        //- Perform the coupling with necessary initialization etc.
        void performCoupling();

        //- Disallow default bitwise copy constructor
        externalCoupled(const externalCoupled&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const externalCoupled&) = delete;


public:

    //- Runtime type information
    TypeName("externalCoupled");

    //- Name of patch key, e.g. '// Patch:' when looking for start of patch data
    static string patchKey;

    //- Inherited variable for logging
    using functionObject::log;


    // Constructors

        //- Construct given time and dictionary
        externalCoupled
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~externalCoupled() = default;


    // Member Functions

      // Function object control

        //- Called at each ++ or += of the time-loop
        virtual bool execute();

        //- Manual execute (sub-loop or when converged)
        virtual bool execute(const label subIndex);

        //- Called when Time::run() determines that the time-loop exits
        virtual bool end();

        //- Read and set the function object if its data have changed
        virtual bool read(const dictionary& dict);

        //- Write, currently a no-op
        virtual bool write();


      // File creation, removal

        //- Write data files (all regions, all fields) from master (OpenFOAM)
        virtual void writeDataMaster() const;

        //- Read data files (all regions, all fields) on master (OpenFOAM)
        virtual void readDataMaster();

        //- Remove data files written by master (OpenFOAM)
        virtual void removeDataMaster() const;

        //- Remove data files written by slave (external code)
        virtual void removeDataSlave() const;


      // Other

        //- Create single name by appending words (in sorted order),
        //  separated by '_'
        static word compositeName(const wordList&);

        //- Write geometry for the group as region/patch
        static void writeGeometry
        (
            const UPtrList<const fvMesh>& meshes,
            const fileName& commsDir,
            const wordRe& groupName
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "externalCoupledTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
