/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "foamChemistryReader.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::foamChemistryReader<ThermoType>::foamChemistryReader
(
    const fileName& reactionsFileName,
    const fileName& thermoFileName
)
:
    chemistryReader<ThermoType>(),
    speciesThermo_(dictionary(IFstream(thermoFileName.expand())())),
    speciesTable_
    (
        dictionary(IFstream(reactionsFileName.expand())()).lookup("species")
    ),
    reactions_
    (
        speciesTable_,
        speciesThermo_,
        dictionary(IFstream(reactionsFileName))
    )
{}


template<class ThermoType>
Foam::foamChemistryReader<ThermoType>::foamChemistryReader
(
    const dictionary& thermoDict
)
:
    chemistryReader<ThermoType>(),
    speciesThermo_
    (
        dictionary
        (
            IFstream
            (
                fileName(thermoDict.lookup("foamChemistryThermoFile")).expand()
            )()
        )
    ),
    speciesTable_
    (
        dictionary
        (
            IFstream
            (
                fileName(thermoDict.lookup("foamChemistryFile")).expand()
            )()
        ).lookup("species")
    ),
    reactions_
    (
        speciesTable_,
        speciesThermo_,
        dictionary
        (
            IFstream
            (
                fileName(thermoDict.lookup("foamChemistryFile")).expand()
            )()
        )
    )
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
