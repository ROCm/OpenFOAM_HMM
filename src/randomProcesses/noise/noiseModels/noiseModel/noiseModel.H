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
    Foam::noiseModel

Description
    Base class for noise models.

    Data is read from a dictionary, e.g.

    \verbatim
    pRef            0;
    N               4096;
    fl              25;
    fu              25;
    startTime       0;
    \endverbatim

    where
    \table
        Property    | Description                   | Required  | Default value
        pRef        | Reference pressure            | no        | 0
        N           | Number of samples in sampling window | no | 65536 (2^16)
        fl          | Lower frequency bounds        | no        | 25
        fu          | Upper frequency bounds        | no        | 10000
    \endtable

Note
    The number of samples in the sampling window must be a power of 2


SourceFiles
    noiseModel.C

\*---------------------------------------------------------------------------*/

#ifndef noiseModel_H
#define noiseModel_H

#include "dictionary.H"
#include "scalarList.H"
#include "instantList.H"
#include "windowModel.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class noiseModel Declaration
\*---------------------------------------------------------------------------*/

class noiseModel
{

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        noiseModel(const noiseModel&);

        //- Disallow default bitwise assignment
        void operator=(const noiseModel&);


protected:

    // Protected Data

        //- Copy of dictionary used for construction
        const dictionary dict_;

        //- Reference pressure
        scalar pRef_;

        //- Number of samples in sampling window, default = 2^16
        label nSamples_;

        //- Lower frequency limit, default = 25Hz
        scalar fLower_;

        //- Upper frequency limit, default = 10kHz
        scalar fUpper_;

        //- Start time, default = 0s
        scalar startTime_;

        //- Window model
        autoPtr<windowModel> windowModelPtr_;
        


    // Protected Member Functions

        //- Check and return uniform time step
        scalar checkUniformTimeStep
        (
            const scalarList& times
        ) const; 

        //- Find and return start time index
        label findStartTimeIndex
        (
            const instantList& allTimes,
            const scalar startTime
        ) const;


public:

    //- Runtime type information
    TypeName("noiseModel");

    //- Run time selection table  
    declareRunTimeSelectionTable
    (
        autoPtr,
        noiseModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    //- Selector
    static autoPtr<noiseModel> New(const dictionary& dict);

    //- Constructor
    noiseModel(const dictionary& dict);

    //- Destructor
    virtual ~noiseModel();


    // Public Member Functions

        //- Abstract call to calculate
        virtual void calculate() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
