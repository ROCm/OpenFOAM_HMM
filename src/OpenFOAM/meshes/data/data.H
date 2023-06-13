/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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
    Foam::data

Description
    Database for solution data, solver performance and other reduced data.

    Both faMesh and fvMesh are derived from this class so that every
    GeometricField has access to this data via their mesh reference.

SourceFiles
    data.C
    dataTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_data_H
#define Foam_data_H

#include "IOdictionary.H"
#include "solverPerformance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class data Declaration
\*---------------------------------------------------------------------------*/

class data
:
    public IOdictionary
{
    // Private Data

        //- Previously used time-index, used for reset between iterations
        mutable label prevTimeIndex_;


    // Private Member Functions

        //- True if boolean entry exists and is set
        bool getBoolEntry(const word& keyword) const
        {
            bool value = false;
            return readIfPresent(keyword, value) && value;
        }

        //- Add/remove boolean entry
        void setBoolEntry(const word& keyword, bool on)
        {
            if (on)
            {
                add(keyword, true);
            }
            else
            {
                remove(keyword);
            }
        }

        //- No copy construct
        data(const data&) = delete;

        //- No copy assignment
        void operator=(const data&) = delete;

public:

    // Constructors

        //- Construct for objectRegistry (registered with specified name)
        //- optionally with initial content
        data
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary* content = nullptr
        );

        //- Construct for objectRegistry (registered as "data")
        explicit data(const objectRegistry& obr)
        :
            data("data", obr)
        {}

        //- Construct for objectRegistry (registered with specified name)
        //- copying initial content
        data
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        )
        :
            data(name, obr, &dict)
        {}

        //- Construct for objectRegistry (registered as "data")
        //- copying initial content
        data(const objectRegistry& obr, const dictionary& dict)
        :
            data("data", obr, &dict)
        {}


    // Member Functions

        //- Return the dictionary of solver performance data which
        //- includes initial and final residuals for convergence checking
        const dictionary& solverPerformanceDict() const;

        //- True if "firstIteration" entry exists and is set
        bool isFirstIteration() const
        {
            return getBoolEntry("firstIteration");
        }

        //- True if "finalIteration" entry exists and is set
        bool isFinalIteration() const
        {
            return getBoolEntry("finalIteration");
        }

        //- Add/remove "firstIteration" entry
        void isFirstIteration(bool on)
        {
            return setBoolEntry("firstIteration", on);
        }

        //- Add/remove "finalIteration" entry
        void isFinalIteration(bool on)
        {
            return setBoolEntry("finalIteration", on);
        }


        //- Add/set the solverPerformance entry for the named field
        template<class Type>
        void setSolverPerformance
        (
            const word& name,
            const SolverPerformance<Type>& sp
        ) const;

        //- Add/set the solverPerformance entry, using its fieldName
        template<class Type>
        void setSolverPerformance
        (
            const SolverPerformance<Type>& sp
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dataTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
