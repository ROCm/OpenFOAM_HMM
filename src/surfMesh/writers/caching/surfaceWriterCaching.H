/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::surfaceWriters::writerCaching

Description
    Information for surface writers with collated times.

    The class maintains an internal list of the known times
    as well as a file-cached version with the field information.
    The information is used for restarts.

SourceFiles
    surfaceWriterCaching.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceWriters_writerCaching_H
#define surfaceWriters_writerCaching_H

#include "bitSet.H"
#include "dictionary.H"
#include "scalarList.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{

/*---------------------------------------------------------------------------*\
                        Class writerCaching Declaration
\*---------------------------------------------------------------------------*/

class writerCaching
{
    // Private Data

        //- Cache dictionary file name
        word dictName_;

        //- The output times
        DynamicList<scalar> times_;

        //- Indices in times_ when geometry (mesh) has been written
        bitSet geoms_;

        //- Cached information for geometry, times, fields
        dictionary cache_;


    // Private Member Functions

        //- Read time information from dictFileName.
        //  Returns timeIndex corresponding to timeValue
        label readPreviousTimes
        (
            const fileName& dictFile,
            const scalar timeValue
        );

        //- Get or create a sub-dictionary for named field
        dictionary& fieldDict(const word& fieldName);

        //- Remove named field
        bool remove(const word& fieldName);


public:

    // Constructors

        //- Construct with specified cache name
        explicit writerCaching(const word& cacheFileName);


    //- Destructor
    virtual ~writerCaching() = default;


    // Member Functions

        //- The output times for fields
        const scalarList& times() const
        {
            return times_;
        }

        //- Indices in times() when geometry (mesh) has been written
        const bitSet& geometries() const
        {
            return geoms_;
        }

        //- The most current time index
        label latestTimeIndex() const
        {
            return max(0, times_.size()-1);
        }

        //- The most current geometry index
        label latestGeomIndex() const
        {
            return max(0, geoms_.find_last());
        }

        //- Get or create the 'fields' information dictionary.
        const dictionary& fieldsDict() const;

        //- Clear all values
        void clear();

        //- Update time/geometry information and file cache.
        //- This routine should only be called from the master process
        //  \return True if there is a state change, which is either a
        //      geometry change or a new time interval
        bool update
        (
            const fileName& baseDir,    //!< Directory containing the cache file
            const scalar timeValue,     //!< The current time value
            const bool geomChanged,     //!< Monitored geometry changed
            const word& fieldName,      //!< Name of field
            const word& fieldType,      //!< Type of field
            const word& varName = word::null //!< Alternative field name
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceWriters
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
