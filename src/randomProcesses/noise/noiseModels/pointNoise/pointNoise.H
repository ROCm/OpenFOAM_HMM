/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenCFD Ltd.
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
    Foam::noiseModels::pointNoise

Description
    Perform noise analysis on point-based pressure data.

    Input data is read from a dictionary, e.g.

    \verbatim
    // Pressure reference
    pRef            0;

    // Number of samples in sampling window
    // Must be a power of 2, default = 2^16 (=65536)
    N               4096;

    // Lower frequency bounds 
    fl              25;

    // Upper frequency bounds 
    fu              25;

    // Start time
    startTime       0;

    windowModel     <modelType>
    <modelType>Coeffs
    {
        ...
    }

    // Pressure data supplied in CSV file format
    csvFileData
    {
        fileName        "pressureData";
        nHeaderLine     1;
        refColumn       0;
        componentColumns (1);
        separator       " ";
        mergeSeparators yes;
    }

    graphFormat     raw;

    \endverbatim

SourceFiles
    pointNoise.C

SeeAlso
    noiseModel.H
    windowModel.H

\*---------------------------------------------------------------------------*/

#ifndef pointNoise_H
#define pointNoise_H

#include "noiseModel.H"
#include "CSV.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace noiseModels
{

/*---------------------------------------------------------------------------*\
                         Class pointNoise Declaration
\*---------------------------------------------------------------------------*/

class pointNoise
:
    public noiseModel
{

protected:

    // Protected data

        //- Graph format
        word graphFormat_;


    // Protected Member Functions

        void filterTimeData
        (
            const CSV<scalar>& pData,
            scalarField& t,
            scalarField& p
        );


public:

    //- Runtime type information
    TypeName("pointNoise");

    //- Constructor
    pointNoise(const dictionary& dict);

    //- Destructor
    virtual ~pointNoise();


    // Public Member Functions

        //- Calculate
        virtual void calculate();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace noiseModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
