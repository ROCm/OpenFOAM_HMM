/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
     \\/     M anipulation  |
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
    fileMonitor

\*----------------------------------------------------------------------------*/

#include "fileMonitor.H"
#include "IOstreams.H"
#include "Pstream.H"
#include "PackedList.H"
#include "PstreamReduceOps.H"

#ifdef FOAM_USE_STAT
#   include "OSspecific.H"
#   include "regIOobject.H"     // for fileModificationSkew symbol
#else
#   include <sys/inotify.h>
#   include <stropts.h>
#   include <sys/ioctl.h>
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::fileMonitor, 0);

template<>
const char* Foam::NamedEnum<Foam::fileMonitor::fileState, 3>::names[] =
{
    "unmodified",
    "deleted",
    "modified"
};
const Foam::NamedEnum<Foam::fileMonitor::fileState, 3>
    Foam::fileMonitor::fileStateNames_;


namespace Foam
{
    // Reduction operator for PackedList of fileState
    class reduceFileStates
    {
        public:
        unsigned int operator()(const unsigned int x, const unsigned int y)
        const
        {
            // x,y are sets of 2bits representing fileState

            unsigned int mask = 3u;
            unsigned int shift = 0;
            unsigned int result = 0;

            while (mask)
            {
                // Combine state
                unsigned int xState = (x & mask) >> shift;
                unsigned int yState = (y & mask) >> shift;

                // Combine and add to result. Combine is such that UNMODIFIED
                // wins.
                unsigned int state = min(xState, yState);
                result |= (state << shift);

                shift += 2;
                mask <<= 2;
            }
            return result;
        }
    };

    // Combine operator for PackedList of fileState
    class combineReduceFileStates
    {
        public:
        void operator()(unsigned int& x, const unsigned int y) const
        {
            x = reduceFileStates()(x, y);
        }
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

#ifdef FOAM_USE_STAT
void Foam::fileMonitor::checkFiles() const
{
    for
    (
        HashTable<label, time_t>::iterator iter = lastModified_.begin();
        iter != lastModified_.end();
        ++iter
    )
    {
        label watchFd = iter.key();
        const fileName& fName = watchFile_[watchFd];
        time_t newTime = lastModified(fName);

        if (newTime == 0)
        {
            state_.set(watchFd, DELETED);
        }
        else
        {
            time_t oldTime = iter();
            if (newTime > (oldTime + regIOobject::fileModificationSkew))
            {
                iter() = newTime;
                state_.set(watchFd, MODIFIED);
            }
            else
            {
                state_.set(watchFd, UNMODIFIED);
            }
        }
    }
}
#else
void Foam::fileMonitor::checkFiles() const
{
    while (true)
    {
        struct timeval zeroTimeout = {0, 0};

        int ready = select
        (
            inotifyFd_+1,       // num filedescriptors in watchSet_
            &watchSet_,         // watchSet_ with only inotifyFd
            NULL,
            NULL,
            &zeroTimeout
        );

        if (ready < 0)
        {
            FatalErrorIn("fileMonitor::updateStates()")
                << "Problem in issuing select."
                << abort(FatalError);
        }
        else if (FD_ISSET(inotifyFd_, &watchSet_))
        {
            struct inotify_event inotifyEvent;

            // Read first event
            ssize_t nBytes = read
            (
                inotifyFd_,
                &inotifyEvent,
                sizeof(inotifyEvent)
            );

            if (nBytes != sizeof(inotifyEvent))
            {
                FatalErrorIn("fileMonitor::updateStates(const fileName&)")
                    << "Read " << label(nBytes) << " ; expected "
                    << label(sizeof(inotifyEvent))
                    << abort(FatalError);
            }

            //Pout<< "mask:" << inotifyEvent.mask << endl;
            //Pout<< "watchFd:" << inotifyEvent.wd << endl;
            //Pout<< "watchName:" << watchFile_[inotifyEvent.wd] << endl;

            switch (inotifyEvent.mask)
            {
                case IN_DELETE_SELF:
                {
                    Map<fileState>::iterator iter =
                        state_.find(label(inotifyEvent.wd));
                    iter() = DELETED;
                }
                break;

                case IN_MODIFY:
                case IN_CLOSE_WRITE:
                {
                    Map<fileState>::iterator iter =
                        state_.find(label(inotifyEvent.wd));
                    iter() = MODIFIED;
                }
                break;
            }
        }
        else
        {
            // No data. Reset watchSet_
            FD_SET(inotifyFd_, &watchSet_);
            return;
        }
    }
}
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
#ifdef FOAM_USE_STAT

Foam::fileMonitor::fileMonitor()
:
    state_(20),
    watchFile_(20),
    lastModified_(20)
{}

#else

Foam::fileMonitor::fileMonitor()
:
    state_(20),
    watchFile_(20),
    inotifyFd_(inotify_init())
{
    //- Add notify descriptor to select set
    FD_ZERO(&watchSet_);
    FD_SET(inotifyFd_, &watchSet_);
}

#endif

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileMonitor::~fileMonitor()
{
    // Remove any remaining files
    List<label> watchFds(state_.toc());
    forAll(watchFds, i)
    {
        removeWatch(watchFds[i]);
    }
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::fileMonitor::addWatch(const fileName& fName)
{
#ifdef FOAM_USE_STAT
    label watchFd = lastModified_.size();
    lastModified_.insert(watchFd, lastModified(fName));
#else
    label watchFd = inotify_add_watch
    (
        inotifyFd_,
        fName.c_str(),
        //IN_ALL_EVENTS
        IN_CLOSE_WRITE | IN_DELETE_SELF | IN_MODIFY
    );
#endif

    if (debug)
    {
        Pout<< "fileMonitor : added watch " << watchFd << " on file "
            << fName << endl;
    }

    if (watchFd < 0)
    {
        WarningIn("fileMonitor::addWatch(const fileName&)")
            << "could not add watch for file " << fName << endl;
    }
    else
    {
        state_.insert(watchFd, UNMODIFIED);
        watchFile_.insert(watchFd, fName);
    }
    return watchFd;
}


bool Foam::fileMonitor::removeWatch(const label watchFd)
{
    if (debug)
    {
        Pout<< "fileMonitor : removing watch " << watchFd << " on file "
            << watchFile_[watchFd] << endl;
    }

    state_.erase(watchFd);
    watchFile_.erase(watchFd);
#ifdef FOAM_USE_STAT
    return lastModified_.erase(watchFd);
#else
    return inotify_rm_watch(inotifyFd_, int(watchFd)) == 0;
#endif
}


const Foam::fileName& Foam::fileMonitor::getFile(const label watchFd) const
{
    return watchFile_[watchFd];
}


Foam::fileMonitor::fileState Foam::fileMonitor::getState(const label watchFd)
const
{
    return state_[watchFd];
}


void Foam::fileMonitor::updateStates(const bool syncPar) const
{
    checkFiles();

    if (syncPar)
    {
        PackedList<2> stats(state_.size());
        label i = 0;
        forAllConstIter(Map<fileState>, state_, iter)
        {
            stats[i++] = (unsigned int)(iter());
        }
        // Save local state for warning message below
        PackedList<2> thisProcStats(stats);

        if (stats.storage().size() == 1)
        {
            // Optimisation valid for most cases.
            reduce(stats.storage()[0], reduceFileStates());
        }
        else
        {
            Pstream::listCombineGather
            (
                stats.storage(),
                combineReduceFileStates()
            );
        }

        i = 0;
        forAllIter(Map<fileState>, state_, iter)
        {
            if (thisProcStats[i] != UNMODIFIED)
            {
                if (stats[i] == UNMODIFIED)
                {
                    WarningIn("fileMonitor::updateStates(const bool) const")
                        << "Delaying reading " << watchFile_[iter.key()]
                        << " due to inconsistent "
                           "file time-stamps between processors"
                        << endl;
                }
                else
                {
                    unsigned int stat = stats[i];
                    iter() = fileState(stat);
                }
            }
            i++;
        }
    }
}


void Foam::fileMonitor::setUnmodified(const label watchFd)
{
#ifdef FOAM_USE_STAT
    lastModified_[watchFd] = lastModified(watchFile_[watchFd]);
#endif

    Map<fileState>::iterator iter = state_.find(watchFd);

    if (iter == state_.end())
    {
        FatalErrorIn("fileMonitor::setUnmodified(const label)")
            << "Illegal watchFd " << watchFd
            << abort(FatalError);
    }

    iter() = UNMODIFIED;
}


// ************************************************************************* //
