/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Application
    mag

Description
    Calculates and writes the magnitude-squared of a field for each time

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void writeMagSqrField
(
    const IOobject& header,
    const fvMesh& mesh,
    bool& processed
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (header.headerClassName() == fieldType::typeName)
    {
        Info<< "    Reading " << header.name() << endl;
        fieldType field(header, mesh);

        Info<< "    Calculating magSqr" << header.name() << endl;
        volScalarField magSqrField
        (
            IOobject
            (
                "magSqr" + header.name(),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ
            ),
            magSqr(field)
        );
        magSqrField.write();

        processed = true;
    }
}


int main(int argc, char *argv[])
{
    argList::validArgs.append("fieldName");

#   include "addTimeOptions.H"
#   include "setRootCase.H"

    word fieldName(args.additionalArgs()[0]);

#   include "createTime.H"

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

#   include "createMesh.H"

    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        IOobject fieldHeader
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check field "fieldName" exists
        if (fieldHeader.headerOk())
        {
            mesh.readUpdate();

            bool processed = false;
            writeMagSqrField<scalar>(fieldHeader, mesh, processed);
            writeMagSqrField<vector>(fieldHeader, mesh, processed);
            writeMagSqrField<sphericalTensor>(fieldHeader, mesh, processed);
            writeMagSqrField<symmTensor>(fieldHeader, mesh, processed);
            writeMagSqrField<tensor>(fieldHeader, mesh, processed);
            if (!processed)
            {
                FatalError
                    << "Unable to process " << fieldName << nl
                    << "No call to mag for fields of type "
                    << fieldHeader.headerClassName() << nl << nl
                    << exit(FatalError);
            }
        }
        else
        {
            Info<< "    No " << fieldName << endl;
        }

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
