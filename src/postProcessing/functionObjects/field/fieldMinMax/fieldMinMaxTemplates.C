/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "fieldMinMax.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::fieldMinMax::calcMinMaxFields
(
    const word& fieldName,
    const modeType& mode
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (obr_.foundObject<fieldType>(fieldName))
    {
        const label procI = Pstream::myProcNo();

        const fieldType& field = obr_.lookupObject<fieldType>(fieldName);
        switch (mode)
        {
            case mdMag:
            {
                const scalarField magField(mag(field));

                labelList minIs(Pstream::nProcs());
                scalarList minVs(Pstream::nProcs());
                List<vector> minCs(Pstream::nProcs());
                minIs[procI] = findMin(magField);
                minVs[procI] = magField[minIs[procI]];
                minCs[procI] = field.mesh().C()[minIs[procI]];

                Pstream::gatherList(minIs);
                Pstream::gatherList(minVs);
                Pstream::gatherList(minCs);

                labelList maxIs(Pstream::nProcs());
                scalarList maxVs(Pstream::nProcs());
                List<vector> maxCs(Pstream::nProcs());
                maxIs[procI] = findMax(magField);
                maxVs[procI] = magField[maxIs[procI]];
                maxCs[procI] = field.mesh().C()[maxIs[procI]];

                Pstream::gatherList(maxIs);
                Pstream::gatherList(maxVs);
                Pstream::gatherList(maxCs);

                if (Pstream::master())
                {
                    label minI = findMin(minVs);
                    scalar minValue = minVs[minI];
                    const vector& minC = minCs[minI];

                    label maxI = findMax(maxVs);
                    scalar maxValue = maxVs[maxI];
                    const vector& maxC = maxCs[maxI];

                    if (write_)
                    {
                        fieldMinMaxFilePtr_()
                            << obr_.time().value() << token::TAB
                            << fieldName << token::TAB
                            << minValue << token::TAB << minC;

                        if (Pstream::parRun())
                        {
                            fieldMinMaxFilePtr_() << token::TAB << minI;
                        }

                        fieldMinMaxFilePtr_()
                             << token::TAB << maxValue << token::TAB << maxC;

                        if (Pstream::parRun())
                        {
                            fieldMinMaxFilePtr_() << token::TAB << maxI;
                        }

                        fieldMinMaxFilePtr_() << endl;
                    }

                    if (log_)
                    {
                        Info<< "    min(mag(" << fieldName << ")) = "
                            << minValue << " at position " << minC;

                        if (Pstream::parRun())
                        {
                            Info<< " on processor " << minI;
                        }

                        Info<< nl << "    max(mag(" << fieldName << ")) = "
                            << maxValue << " at position " << maxC;

                        if (Pstream::parRun())
                        {
                            Info<< " on processor " << maxI;
                        }

                        Info<< endl;
                    }
                }
                break;
            }
            case mdCmpt:
            {
                List<Type> minVs(Pstream::nProcs());
                labelList minIs(Pstream::nProcs());
                List<vector> minCs(Pstream::nProcs());
                minIs[procI] = findMin(field);
                minVs[procI] = field[minIs[procI]];
                minCs[procI] = field.mesh().C()[minIs[procI]];

                Pstream::gatherList(minIs);
                Pstream::gatherList(minVs);
                Pstream::gatherList(minCs);

                List<Type> maxVs(Pstream::nProcs());
                labelList maxIs(Pstream::nProcs());
                List<vector> maxCs(Pstream::nProcs());
                maxIs[procI] = findMax(field);
                maxVs[procI] = field[maxIs[procI]];
                maxCs[procI] = field.mesh().C()[maxIs[procI]];

                Pstream::gatherList(maxIs);
                Pstream::gatherList(maxVs);
                Pstream::gatherList(maxCs);

                if (Pstream::master())
                {
                    label minI = findMin(minVs);
                    Type minValue = minVs[minI];
                    const vector& minC = minCs[minI];

                    label maxI = findMax(maxVs);
                    Type maxValue = maxVs[maxI];
                    const vector& maxC = maxCs[maxI];

                    if (write_)
                    {
                        fieldMinMaxFilePtr_()
                            << obr_.time().value() << token::TAB
                            << fieldName << token::TAB
                            << minValue << token::TAB << minC;

                        if (Pstream::parRun())
                        {
                            fieldMinMaxFilePtr_() << token::TAB << minI;
                        }

                        fieldMinMaxFilePtr_()
                             << token::TAB << maxValue << token::TAB << maxC;

                        if (Pstream::parRun())
                        {
                            fieldMinMaxFilePtr_() << token::TAB << maxI;
                        }

                        fieldMinMaxFilePtr_() << endl;
                    }

                    if (log_)
                    {
                        Info<< "    min(" << fieldName << ") = "
                            << minValue << " at position " << minC;

                        if (Pstream::parRun())
                        {
                            Info<< " on processor " << minI;
                        }

                        Info<< nl << "    max(" << fieldName << ") = "
                            << maxValue << " at position " << maxC;

                        if (Pstream::parRun())
                        {
                            Info<< " on processor " << maxI;
                        }

                        Info<< endl;
                    }
                }
                break;
            }
            default:
            {
                FatalErrorIn("Foam::fieldMinMax::calcMinMaxFields(const word&)")
                    << "Unknown min/max mode: " << modeTypeNames_[mode_]
                    << exit(FatalError);
            }
        }
    }
}


// ************************************************************************* //
