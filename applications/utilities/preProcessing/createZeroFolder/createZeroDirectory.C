/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

Application
    createZeroDirectory

Description
    Creates a zero directory with fields appropriate for the chosen solver and
    turbulence model. Operates on both single and multi-region cases.

Usage
    The set-up is configured using a 'caseProperties' dictionary, located under
    the $FOAM_CASE/system (or system/regionName if multi-region) directory.
    This consists of a lists of initial and boundary conditions, e.g.

    \verbatim
    initialConditions
    {
        U           uniform (0 0 0);
        p           uniform 0;
    }

    boundaryConditions
    {
        topWall
        {
            category        wall;
            patches         (movingWall);
            type            noSlip;
            options
            {
                wallFunction    highReynolds;
                motion          moving;
            };
            values
            {
                U           uniform (1 0 0);
            }
        }

        walls
        {
            category        wall;
            patches         (fixedWalls);
            type            noSlip;
            options
            {
                wallFunction    highReynolds;
                motion          stationary;
            };
        }
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "boundaryInfo.H"
#include "caseInfo.H"
#include "boundaryTemplates.H"
#include "solverTemplate.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const word getClassType(const word& pType)
{
    if (pType == pTraits<scalar>::typeName)
    {
        return GeometricField<scalar, fvPatchField, volMesh>::typeName;
    }
    else if (pType == pTraits<vector>::typeName)
    {
        return GeometricField<vector, fvPatchField, volMesh>::typeName;
    }
    else if (pType == pTraits<sphericalTensor>::typeName)
    {
        return GeometricField<sphericalTensor, fvPatchField, volMesh>::typeName;
    }
    else if (pType == pTraits<symmTensor>::typeName)
    {
        return GeometricField<symmTensor, fvPatchField, volMesh>::typeName;
    }
    else if (pType == pTraits<tensor>::typeName)
    {
        return GeometricField<tensor, fvPatchField, volMesh>::typeName;
    }
    else
    {
        // error
        return word::null;
    }
}


void createFieldFiles
(
    const Time& runTime,
    const word& regionName,
    const word& regionPrefix,
    const wordList& fieldNames,
    const wordList& fieldTypes,
    const PtrList<dimensionSet>& fieldDimensions,
    const caseInfo& cInfo,
    const boundaryTemplates& bcTemplates
)
{
    Info<< "    Generating field files" << nl << endl;

    // create files
    mkDir(runTime.path()/"0"/regionName);
    forAll(fieldNames, i)
    {
        const_cast<word&>(IOdictionary::typeName) =
            getClassType(fieldTypes[i]);

//        IOdictionary field
//        (
//            IOobject
//            (
//                fieldNames[i],
//                "0",
//                regionName,
//                runTime,
//                IOobject::NO_READ
//            )
//        );

        dictionary field;

        word regionPath = "/";

        if (regionName != word::null)
        {
            regionPath += regionName + '/';
        }

        field.add
        (
            "#include",
            string
            (
                "${FOAM_CASE}/system" + regionPath + "caseProperties"
            )
        );

        field.add("dimensions", fieldDimensions[i]);

        string iField("${:initialConditions." + fieldNames[i] + '}');
        field.add("internalField", iField.c_str());

        dictionary boundaryField =
            cInfo.generateBoundaryField
            (
                regionPrefix,
                fieldNames[i],
                bcTemplates
            );

        field.add("boundaryField", boundaryField);

        // expand all of the dictionary redirections and remove unnecessary
        // entries
        OStringStream os;
        os  << field;

        entry::disableFunctionEntries = 0;
        dictionary field2(IStringStream(os.str())());
        entry::disableFunctionEntries = 1;
        field2.remove("#include");
        field2.remove("initialConditions");
        field2.remove("boundaryConditions");

        // construct and write field dictionary
        IOdictionary fieldOut
        (
            IOobject
            (
                fieldNames[i],
                "0",
                regionName,
                runTime,
                IOobject::NO_READ
            ),
            field2
        );

        fieldOut.regIOobject::writeObject
        (
            IOstream::ASCII,
            IOstream::currentVersion,
            IOstream::UNCOMPRESSED
        );
    }
}


// Main program:
int main(int argc, char *argv[])
{
    // keep variable substitutions
    entry::disableFunctionEntries = 1;

    #include "setRootCase.H"
    #include "createTime.H"

    Info<< "Reading controlDict" << nl << endl;

    IOdictionary controlDict
    (
        IOobject
        (
            "controlDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    fileName baseDir("${WM_PROJECT_DIR}/etc/templates");
    baseDir.expand();

    // read the solver
    const word& solverName = controlDict.lookup("application");

    // generate solver template
    const solverTemplate solver(baseDir, runTime, solverName);

    // read the boundary condition templates
    const boundaryTemplates bcTemplates
    (
        baseDir,
        runTime,
        solver.type()
    );

    Info<< endl;

    const label nRegion = solver.nRegion();
    for (label regionI = 0; regionI < nRegion; regionI++)
    {
        const word& regionName = solver.regionName(regionI);

        if (regionName == word::null)
        {
            Info<< "Region: " << polyMesh::defaultRegion << " (default)"
                << endl;
        }
        else
        {
            Info<< "Region: " << regionName << endl;
        }

        // read the case set-up info for the current region
        const caseInfo cInfo(runTime, regionName);

        cInfo.checkPatches(solver.regionType(regionI), bcTemplates);

        createFieldFiles
        (
            runTime,
            regionName,
            solver.regionType(regionI),
            solver.fieldNames(regionI),
            solver.fieldTypes(regionI),
            solver.fieldDimensions(regionI),
            cInfo,
            bcTemplates
        );
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
