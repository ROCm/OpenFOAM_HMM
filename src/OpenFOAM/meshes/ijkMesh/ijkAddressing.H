/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::ijkAddressing

Description
    A simple i-j-k (row-major order) to linear addressing.

SourceFiles
    ijkAddressingI.H

\*---------------------------------------------------------------------------*/

#ifndef ijkAddressing_H
#define ijkAddressing_H

#include "labelVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ijkAddressing Declaration
\*---------------------------------------------------------------------------*/

class ijkAddressing
{
    // Private Data

        //- The number of cells in the i,j,k directions.
        labelVector sizes_;


public:

    // Constructors

        //- Construct zero-size addressing
        inline ijkAddressing();

        //- Construct with addressing
        inline explicit ijkAddressing(const labelVector& ijk);

        //- Construct with addressing components
        inline ijkAddressing(const label ni, const label nj, const label nk);


    // Member Functions

    // Access

        //- Return the i,j,k addressing sizes
        inline const labelVector& sizes() const;

        //- Return the i,j,k addressing sizes for modification
        inline labelVector& sizes();

        //- Return the total i*j*k size
        inline label size() const;

        //- Addressing is considered empty if any component is zero
        inline bool empty() const;

        //- Reset to (0,0,0) sizing
        inline void clear();

        //- Change the sizing parameters
        inline void reset(const label ni, const label nj, const label nk);

        //- Linear addressing index (offset) for an i,j,k position.
        inline label index(const label i, const label j, const label k) const;


    // Checks

        //- Check indices are within valid ni,nj,nk range.
        //  Optionally allow an extra index for point addressing
        inline void checkIndex
        (
            const label i,
            const label j,
            const label k,
            const bool isPoint = false
        ) const;

        //- Check that sizes() are valid
        inline void checkSizes() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ijkAddressingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
