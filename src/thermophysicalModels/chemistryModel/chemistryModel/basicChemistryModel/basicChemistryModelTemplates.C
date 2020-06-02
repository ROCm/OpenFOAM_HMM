/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "basicChemistryModel.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::autoPtr<ChemistryModel> Foam::basicChemistryModel::New
(
    typename ChemistryModel::reactionThermo& thermo
)
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
            false // Do not register
        )
    );

    if (!chemistryDict.isDict("chemistryType"))
    {
        FatalErrorInFunction
            << "Template parameter based chemistry solver selection is no "
            << "longer supported. Please create a chemistryType dictionary"
            << "instead." << endl << endl << "For example, the entry:" << endl
            << "    chemistrySolver ode<StandardChemistryModel<"
            << "rhoChemistryModel,sutherlandspecie<janaf<perfectGas>,"
            << "sensibleInternalEnergy>>" << endl << endl << "becomes:" << endl
            << "    chemistryType" << endl << "    {" << endl
            << "        solver ode;" << endl << "        method standard;"
            << endl << "    }" << exit(FatalError);
    }

    const dictionary& chemistryTypeDict =
        chemistryDict.subDict("chemistryType");

    const word solverName
    (
        chemistryTypeDict.getCompat<word>
        (
            "solver",
            {{"chemistrySolver", -1712}}
        )
    );

    const word methodName
    (
        chemistryTypeDict.getOrDefault<word>
        (
            "method",
            chemistryTypeDict.getOrDefault<bool>("TDAC", false)
          ? "TDAC"
          : "standard"
        )
    );

    dictionary chemistryTypeDictNew;
    chemistryTypeDictNew.add("solver", solverName);
    chemistryTypeDictNew.add("method", methodName);

    Info<< "Selecting chemistry solver " << chemistryTypeDictNew << endl;

    auto* cstrTable = ChemistryModel::thermoConstructorTablePtr_;

    const word chemSolverCompThermoName
    (
        solverName + '<' + methodName + '<'
      + ChemistryModel::reactionThermo::typeName + ','
      + thermo.thermoName() + ">>"
    );

    auto cstrIter = cstrTable->cfind(chemSolverCompThermoName);

    if (!cstrIter.found())
    {
        constexpr const int nCmpt = 8;

        wordList thisCmpts;
        thisCmpts.append(word::null);
        thisCmpts.append(word::null);
        thisCmpts.append(ChemistryModel::reactionThermo::typeName);
        thisCmpts.append
        (
            basicThermo::splitThermoName(thermo.thermoName(), 5)
        );

        List<wordList> validNames;

        validNames.append
        (
            // Header
            wordList({"solver", "method"})
        );

        List<wordList> validCmpts;
        validCmpts.append
        (
            // Header
            wordList
            ({
                "solver",
                "method",
                "reactionThermo",
                "transport",
                "thermo",
                "equationOfState",
                "specie",
                "energy"
            })
        );

        for (const word& validName : cstrTable->sortedToc())
        {
            validCmpts.append
            (
                basicThermo::splitThermoName(validName, nCmpt)
            );

            const wordList& cmpts = validCmpts.last();

            bool isValid = true;
            for (label i = 2; i < cmpts.size() && isValid; ++i)
            {
                isValid = isValid && cmpts[i] == thisCmpts[i];
            }

            if (isValid)
            {
                validNames.append(SubList<word>(cmpts, 2));
            }
        }


        FatalErrorInFunction
            << "Unknown " << typeName_() << " type " << solverName << nl << nl;

        FatalErrorInFunction
            << "All " << validNames[0][0] << '/' << validNames[0][1]
            << "combinations for this thermodynamic model:" << nl << nl;

        // Table of available packages (as constituent parts)
        printTable(validNames, FatalErrorInFunction)
            << nl
            << "All " << validCmpts[0][0] << '/' << validCmpts[0][1] << '/'
            << validCmpts[0][2] << "/thermoPhysics combinations are:"
            << nl << nl;

        // Table of available packages (as constituent parts)
        printTable(validCmpts, FatalErrorInFunction)
            << exit(FatalError);
    }

    return autoPtr<ChemistryModel>(cstrIter()(thermo));
}

// ************************************************************************* //
