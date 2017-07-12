/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "chemistryReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::autoPtr<Foam::chemistryReader<ThermoType>>
Foam::chemistryReader<ThermoType>::New
(
    const dictionary& thermoDict,
    speciesTable& species
)
{
    // Let the chemistry reader type default to CHEMKIN
    // for backward compatibility
    word readerName("chemkinReader");

    // otherwise use the specified reader
    thermoDict.readIfPresent("chemistryReader", readerName);

    Info<< "Selecting chemistryReader " << readerName << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(readerName);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown chemistryReader type "
            << readerName << nl << nl
            << "Valid chemistryReader types :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<chemistryReader<ThermoType>>
    (
        cstrIter()(thermoDict, species)
    );
}


// ************************************************************************* //
