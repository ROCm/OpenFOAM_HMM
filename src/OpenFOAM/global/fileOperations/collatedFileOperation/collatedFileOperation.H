/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::fileOperations::collatedFileOperation

Description
    Version of masterUncollatedFileOperation that collates regIOobjects
    into a container in the processors/ subdirectory.

    Uses threading if maxThreadFileBufferSize != 0.
        > 0 : Can use mpi inside thread to collect data if buffer is not
              large enough. Does need full thread support inside MPI.

        < 0 : special : -maxThreadFileBufferSize is guaranteed large enough
              for all writing. Initialises MPI without thread support.

See also
    masterUncollatedFileOperation

SourceFiles
    collatedFileOperation.C

\*---------------------------------------------------------------------------*/

#ifndef fileOperations_collatedFileOperation_H
#define fileOperations_collatedFileOperation_H

#include "masterUncollatedFileOperation.H"
#include "OFstreamCollator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileOperations
{

/*---------------------------------------------------------------------------*\
                    Class collatedFileOperation Declaration
\*---------------------------------------------------------------------------*/

class collatedFileOperation
:
    public masterUncollatedFileOperation
{
protected:

    // Protected Data

        //- Any communicator allocated by me
        const label myComm_;

        //- Threaded writer
        mutable OFstreamCollator writer_;

        // For non-parallel operation

            //- Number of processors (overall)
            label nProcs_;

            //- Ranks of IO handlers
            const labelList ioRanks_;


    // Protected Member Functions

        //- Print banner information, optionally with io ranks
        void printBanner(const bool printRanks = false) const;

        //- Is proci master of communicator (in parallel) or master of
        //- the io ranks (non-parallel)
        bool isMasterRank(const label proci) const;

        //- Append to processorsNN/ file
        bool appendObject
        (
            const regIOobject& io,
            const fileName& pathName,
            IOstreamOption streamOpt
        ) const;


public:

        //- Runtime type information
        TypeName("collated");


    // Static Data

        //- Max size of thread buffer size. This is the overall size of
        //  all files to be written. Starts blocking if not enough size.
        //  Read as float to enable easy specification of large sizes.
        static float maxThreadFileBufferSize;


    // Constructors

        //- Default construct
        explicit collatedFileOperation(bool verbose);

        //- Construct from user communicator
        collatedFileOperation
        (
            const label comm,
            const labelList& ioRanks,
            const word& typeName,
            bool verbose
        );


    //- Destructor
    virtual ~collatedFileOperation();


    // Member Functions

        // (reg)IOobject functionality

            //- Generate disk file name for object. Opposite of filePath.
            virtual fileName objectPath
            (
                const IOobject& io,
                const word& typeName
            ) const;

            //- Writes a regIOobject (so header, contents and divider).
            //  Returns success state.
            virtual bool writeObject
            (
                const regIOobject&,
                IOstreamOption streamOpt = IOstreamOption(),
                const bool valid = true
            ) const;

        // Other

            //- Forcibly wait until all output done. Flush any cached data
            virtual void flush() const;

            //- Actual name of processors dir
            virtual word processorsDir(const IOobject&) const;

            //- Actual name of processors dir
            virtual word processorsDir(const fileName&) const;

            //- Set number of processor directories/results.
            //- Only used in decomposePar
            virtual void setNProcs(const label nProcs);
};


/*---------------------------------------------------------------------------*\
               Class collatedFileOperationInitialise Declaration
\*---------------------------------------------------------------------------*/

class collatedFileOperationInitialise
:
    public masterUncollatedFileOperationInitialise
{
public:

    // Constructors

        //- Construct from components
        collatedFileOperationInitialise(int& argc, char**& argv)
        :
            masterUncollatedFileOperationInitialise(argc, argv)
        {}


    //- Destructor
    virtual ~collatedFileOperationInitialise() = default;


    // Member Functions

        //- Requires threading for non-zero maxThreadFileBufferSize
        virtual bool needsThreading() const
        {
            return (collatedFileOperation::maxThreadFileBufferSize > 0);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileOperations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
