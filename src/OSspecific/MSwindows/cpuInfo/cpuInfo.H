/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::cpuInfo

Description
    General CPU characteristics.

    If the machine has multiple cpus/cores, only the characteristics
    of the first core are used.

Note
    Windows variant only provides the number of cores.

SourceFiles
    cpuInfo.C

\*---------------------------------------------------------------------------*/

#ifndef cpuInfo_H
#define cpuInfo_H

#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class Ostream;

/*---------------------------------------------------------------------------*\
                           Class cpuInfo Declaration
\*---------------------------------------------------------------------------*/

class cpuInfo
{
    // Private data

        // Various bits from /proc/cpuinfo

        std::string vendor_id;
        std::string model_name;
        int cpu_family;
        int model;
        float cpu_MHz;
        int siblings;
        int cpu_cores;


    // Private Member Functions

        //- Parse /proc/cpuinfo
        void parse();

        //- No copy construct
        cpuInfo(const cpuInfo&) = delete;

        //- No copy assignment
        void operator=(const cpuInfo&) = delete;

public:

    // Constructors

        //- Construct and populate with information
        cpuInfo();

    //- Destructor
    ~cpuInfo() = default;


    // Member Functions

        //- Write content as dictionary entries
        void write(Ostream& os) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
