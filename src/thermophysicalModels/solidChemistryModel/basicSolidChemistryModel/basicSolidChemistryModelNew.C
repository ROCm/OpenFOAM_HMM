/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "basicSolidChemistryModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicSolidChemistryModel>
Foam::basicSolidChemistryModel::New(solidReactionThermo& thermo)
{
    const IOdictionary chemistryDict
    (
        IOobject
        (
            thermo.phasePropertyName("chemistryProperties"),
            thermo.db().time().constant(),
            thermo.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const dictionary& chemistryTypeDict =
        chemistryDict.subDict("chemistryType");

    Info<< "Selecting chemistry type " << chemistryTypeDict << endl;

    std::initializer_list<const char*> cmptNames
    {
        "chemistrySolver",
        "chemistryThermo",
        "baseChemistry",
        "transport",
        "thermo",
        "equationOfState",
        "specie",
        "energy",
        "transport",
        "thermo",
        "equationOfState",
        "specie",
        "energy"
    };

    const IOdictionary thermoDict
    (
        IOobject
        (
            basicThermo::dictName,
            thermo.db().time().constant(),
            thermo.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    const dictionary& solidThermoTypeDict = thermoDict.subDict("thermoType");
    const word solidThermoTypeName
    (
        solidThermoTypeDict.get<word>("transport") + '<'
      + solidThermoTypeDict.get<word>("thermo") + '<'
      + solidThermoTypeDict.get<word>("equationOfState") + '<'
      + solidThermoTypeDict.get<word>("specie") + ">>,"
      + solidThermoTypeDict.get<word>("energy") + ">"
    );

    const dictionary& gasThermoTypeDict = thermoDict.subDict("gasThermoType");
    const word gasThermoTypeName
    (
        gasThermoTypeDict.get<word>("transport") + '<'
      + gasThermoTypeDict.get<word>("thermo") + '<'
      + gasThermoTypeDict.get<word>("equationOfState") + '<'
      + gasThermoTypeDict.get<word>("specie") + ">>,"
      + gasThermoTypeDict.get<word>("energy") + ">"
    );

    // Construct the name of the chemistry type from the components
    const word chemistryTypeName
    (
        chemistryTypeDict.get<word>("chemistrySolver") + '<'
      + chemistryTypeDict.get<word>("chemistryThermo") + '<'
      + typeName + ','
      + solidThermoTypeName + ',' + gasThermoTypeName + ">>"
    );

    Info<< "chemistryTypeName " << chemistryTypeName << endl;

    auto cstrIter = thermoConstructorTablePtr_->cfind(chemistryTypeName);

    if (!cstrIter.found())
    {
        const int nCmpt = cmptNames.size();

        // Build a table of the thermo packages constituent parts
        // Note: row-0 contains the names of constituent parts
        List<wordList> validCmpts(thermoConstructorTablePtr_->size()+1);

        // Header (row 0)
        validCmpts[0].resize(nCmpt);
        std::copy(cmptNames.begin(), cmptNames.end(), validCmpts[0].begin());

        label rowi = 1;
        for (const word& validName : thermoConstructorTablePtr_->sortedToc())
        {
            validCmpts[rowi] = basicThermo::splitThermoName(validName, nCmpt);

            if (validCmpts[rowi].size())
            {
                ++rowi;
            }
        }
        validCmpts.resize(rowi);


        FatalIOErrorInLookup
        (
            chemistryTypeDict,
            typeName,
            word::null, // Suppress long name? Just output dictionary (above)
            *thermoConstructorTablePtr_
        );

        // Table of available packages (as constituent parts)
        printTable(validCmpts, FatalIOError)
            << exit(FatalIOError);
    }

    return
        autoPtr<basicSolidChemistryModel>(cstrIter()(thermo));
}


// ************************************************************************* //
