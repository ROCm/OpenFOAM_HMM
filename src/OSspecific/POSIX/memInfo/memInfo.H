/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::memInfo

Description
    Memory usage information for the process running this object.

Note
    Uses the information from /proc/PID/status

SourceFiles
    memInfo.C

\*---------------------------------------------------------------------------*/

#ifndef memInfo_H
#define memInfo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class memInfo;
class Istream;
class Ostream;

Istream& operator>>(Istream& is, memInfo& m);
Ostream& operator<<(Ostream& os, const memInfo& m);


/*---------------------------------------------------------------------------*\
                           Class memInfo Declaration
\*---------------------------------------------------------------------------*/

class memInfo
{
    // Private data

        //- Peak memory used by the process (VmPeak in /proc/PID/status)
        int peak_;

        //- Memory used by the process (VmSize in /proc/PID/status)
        int size_;

        //- Resident set size of the process (VmRSS in /proc/PID/status)
        int rss_;

public:

    // Constructors

        //- Construct null
        memInfo();


    //- Destructor
    ~memInfo();


    // Member Functions

        //- Update according to /proc/PID/status contents
        const memInfo& update();


      // Access

        //- Peak memory (VmPeak in /proc/PID/status) at last update()
        inline int peak() const
        {
            return peak_;
        }

        //- Memory size (VmSize in /proc/PID/status) at last update()
        inline int size() const
        {
            return size_;
        }

        //- Resident set size (VmRSS in /proc/PID/status) at last update()
        inline int rss() const
        {
            return rss_;
        }

        //- True if the memory information appears valid
        bool valid() const;


      // Write

        //- Write content as dictionary entries
        void write(Ostream& os) const;


    // IOstream Operators

        //- Read peak/size/rss from stream
        friend Istream& operator>>(Istream& is, memInfo& m);

        //- Write peak/size/rss to stream
        friend Ostream& operator<<(Ostream& os, const memInfo& m);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
