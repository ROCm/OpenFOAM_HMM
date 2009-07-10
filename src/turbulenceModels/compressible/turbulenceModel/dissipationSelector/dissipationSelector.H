/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::compressible::dissipationSelector

Description
    Simplistic runtime selectable handling of thermal dissipation

    Example use:
    @code
        dissipationSelector thermalDissipation(turbulence, dictionary);

        fvScalarMatrix hEqn ( ... );

        if (thermalDissipation.enabled())
        {
            hEqn += thermalDissipation.dissipation();
        }
    @endcode

SourceFiles
    dissipationSelector.C

\*---------------------------------------------------------------------------*/

#ifndef compressibledissipationSelector_H
#define compressibledissipationSelector_H

#include "NamedEnum.H"
#include "dictionary.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace compressible
{

/*---------------------------------------------------------------------------*\
                   Class dissipationSelector Declaration
\*---------------------------------------------------------------------------*/

class dissipationSelector
{

public:

    //- Type of treatment
    enum treatment
    {
        none,
        cascade,
        equilibrium
    };

    static const NamedEnum<treatment, 3> typeNames;

protected:

    // Protected data

        const turbulenceModel& turbModel_;

        //- Type of treatment
        treatment treatment_;

public:

    // Declare name of the class and its debug switch
    ClassName("dissipationSelector");


    // Constructors

        //- Construct given treatment type
        dissipationSelector
        (
            const turbulenceModel&,
            const treatment = none
        );

        //- Construct given treatment type
        dissipationSelector
        (
            const turbulenceModel&,
            const dictionary&
        );


    // Member functions

        //- Get the treatment type
        treatment treatmentType() const
        {
            return treatment_;
        }

        //- Set the treatment type, return old value
        treatment treatmentType(const treatment val)
        {
            treatment old = treatment_;
            treatment_ = val;
            return old;
        }


        //- Set the treatment type, return old value
        const char* name() const
        {
            return typeNames[treatment_];
        }

        //- Set the treatment type, return old value
        bool enabled() const
        {
            return treatment_ != none;
        }

        tmp<volScalarField> dissipation() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
