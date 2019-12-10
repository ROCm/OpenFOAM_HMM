/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Application
    setExprFields

Group
    grpPreProcessingUtilities

Description
    Set values on a selected set of cells/patch-faces via a dictionary.

Note
    Based on funkySetFields
    Copyright 2006-2018 Bernhard Gschaider <bgschaid@hfd-research.com>

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "pointMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "exprOps.H"
#include "volumeExprDriver.H"
#include "timeSelector.H"
#include "dlLibraryTable.H"


using namespace Foam;

using FieldAssociation = expressions::volumeExpr::FieldAssociation;

word fieldGeoType(const FieldAssociation geoType)
{
    switch (geoType)
    {
        case FieldAssociation::VOLUME_DATA : return "cells"; break;
        case FieldAssociation::SURFACE_DATA : return "faces"; break;
        case FieldAssociation::POINT_DATA : return "points"; break;
        default: break;
    }

    return "unknown";
}


//- Simple control structure to with collected switches to simplify passing
struct setExprFieldsControl
{
    bool dryRun;
    bool debugParsing;
    bool cacheVariables;
    bool useDimension;
    bool createNew;
    bool keepPatches;
    bool correctPatches;
    bool correctBCs;
};


template<class Type>
void doCorrectBoundaryConditions
(
    bool correctBCs,
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    if (correctBCs)
    {
        Info<< "Correcting boundary conditions: " << field.name() << nl;
        field.correctBoundaryConditions();
    }
}


template<class Type>
void doCorrectBoundaryConditions
(
    bool correctBCs,
    GeometricField<Type, pointPatchField, pointMesh>& field
)
{
    if (correctBCs)
    {
        Info<< "Correcting boundary conditions: " << field.name() << nl;
        field.correctBoundaryConditions();
    }
}


template<class Type>
void doCorrectBoundaryConditions
(
    bool correctBCs,
    GeometricField<Type, fvsPatchField, surfaceMesh>& field
)
{}


template<class GeoField, class Mesh>
void setField
(
    const word& fieldName,
    const Mesh& mesh,
    const GeoField& result,
    const scalarField& cond,
    const dimensionSet& dims,
    const wordList& valuePatches,

    const setExprFieldsControl& ctrl
)
{
    Info<< "setField(" << fieldName << "): "
        << pTraits<GeoField>::typeName << endl;

    tmp<GeoField> toutput;

    if (ctrl.createNew)
    {
        // Create with zero
        toutput = GeoField::New
        (
            fieldName,
            mesh,
            dimensioned<typename GeoField::value_type>(dims)
        );
    }
    else
    {
        // Read
        toutput = tmp<GeoField>::New
        (
            IOobject
            (
                fieldName,
                mesh.thisDb().time().timeName(),
                mesh.thisDb(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false  // No register
            ),
            mesh
        );
    }

    auto& output = toutput.ref();

    label setCells = 0;

    if (cond.empty())
    {
        // No condition
        output = result;

        setCells = output.size();
    }
    else
    {
        forAll(output, celli)
        {
            if (expressions::boolOp<scalar>()(cond[celli]))
            {
                output[celli] = result[celli];
                ++setCells;
            }
        }
    }

    const label totalCells = returnReduce(output.size(), plusOp<label>());
    reduce(setCells, plusOp<label>());

    forAll(result.boundaryField(), patchi)
    {
        auto& pf = output.boundaryFieldRef()[patchi];

        if (pf.patch().coupled())
        {
            pf == result.boundaryField()[patchi];
        }
    }


    if (setCells == totalCells)
    {
        Info<< "Set all ";
    }
    else
    {
        Info<< "Set " << setCells << " of ";
    }
    Info<< totalCells << " cells" << endl;


    doCorrectBoundaryConditions(ctrl.correctBCs, output);

    if (ctrl.useDimension)
    {
        Info<< "Setting dimension to " << dims << endl;
        output.dimensions().reset(dims);
    }

    if (ctrl.dryRun)
    {
        Info<< "(dry-run): Writing to " << output.name() << nl;
    }
    else
    {
        Info<< "Writing to " << output.name() << nl;
        output.write();
    }
}


void evaluate
(
    const fvMesh& mesh,
    const word& fieldName,
    const expressions::exprString& expression,
    const expressions::exprString& condition,
    const dictionary& dict,
    const dimensionSet& dims,
    const wordList& valuePatches,

    const setExprFieldsControl& ctrl
)
{
    word oldFieldType;

    if (ctrl.createNew)
    {
        Info<< "Set new field: " << fieldName;
    }
    else
    {
        IOobject io
        (
            fieldName,
            mesh.thisDb().time().timeName(),
            mesh.thisDb(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
        io.typeHeaderOk<IOobject>(false);

        oldFieldType = io.headerClassName();

        if (oldFieldType == IOobject::typeName)
        {
            FatalErrorInFunction
                << "Field " << fieldName << " is  "
                << oldFieldType
                << ". Seems that it does not exist. Use 'create'"
                << nl
                << exit(FatalError);
        }

        Info<< "Modify field: " << fieldName
            << " (type " << oldFieldType << ')';
    }

    Info<< " time=" << mesh.thisDb().time().timeName() << nl
        << "Expression:" << nl
        << ">>>>" << nl
        << expression.c_str() << nl
        << "<<<<" << nl
        << "Condition:" << nl
        << ">>>>" << nl
        << condition.c_str() << nl
        << "<<<<" << nl;

    if (ctrl.keepPatches)
    {
        Info<< "Keeping patches unaltered" << endl;
    }
    else if (!valuePatches.empty())
    {
        Info<< "Setting patches " << flatOutput(valuePatches)
            << " to fixed value" << endl;
    }

    Info<< endl;

    expressions::volumeExprDriver driver
    (
        mesh,
        ctrl.cacheVariables
    );

    driver.readDict(dict);

    if (ctrl.debugParsing)
    {
        Info<< "Parsing expression: " << expression << "\nand condition "
            << condition << nl << endl;
        driver.setDebugging(true, true);
    }


    driver.clearVariables();

    scalarField conditionField;

    bool evaluatedCondition = false;

    FieldAssociation conditionDataType(FieldAssociation::VOLUME_DATA);

    if (condition.size() && condition != "true")
    {
        if (ctrl.debugParsing)
        {
            Info<< "Parsing condition:" << condition << endl;
        }

        driver.parse(condition);
        if (ctrl.debugParsing)
        {
            Info<< "Parsed condition" << endl;
        }

        // Process any/all scalar fields. May help with diagnosis

        bool goodCond = true;
        while (goodCond)
        {
            // volScalarField
            {
                const auto* ptr = driver.isResultType<volScalarField>();
                if (ptr)
                {
                    conditionField = ptr->internalField();
                    break;
                }
            }

            // surfaceScalarField
            {
                const auto* ptr = driver.isResultType<surfaceScalarField>();
                if (ptr)
                {
                    conditionField = ptr->internalField();
                    conditionDataType = FieldAssociation::SURFACE_DATA;
                    break;
                }
            }

            // pointScalarField
            {
                const auto* ptr = driver.isResultType<pointScalarField>();
                if (ptr)
                {
                    conditionField = ptr->internalField();
                    conditionDataType = FieldAssociation::POINT_DATA;
                    break;
                }
            }

            // No matching field types
            goodCond = false;
        }

        // Verify that it also logical
        goodCond = goodCond && driver.isLogical();

        if (!goodCond)
        {
            FatalErrorInFunction
                << " condition: " << condition
                << " does not evaluate to a logical expression: "
                << driver.resultType() << nl
                #ifdef FULLDEBUG
                << "contents: " << conditionField
                #endif
                << exit(FatalError);
        }

        if (ctrl.debugParsing)
        {
            Info<< "Condition evaluates to "
                << conditionField << nl;
        }

        evaluatedCondition = true;
    }

    if (ctrl.debugParsing)
    {
        Info<< "Parsing expression:" << expression << endl;
    }

    driver.parse(expression);

    if (ctrl.debugParsing)
    {
        Info<< "Parsed expression" << endl;
    }

    if (evaluatedCondition)
    {
        if (conditionDataType != driver.fieldAssociation())
        {
            FatalErrorInFunction
                << "Mismatch between condition geometric type ("
                << fieldGeoType(conditionDataType) << ") and" << nl
                << "expression geometric type ("
                << fieldGeoType(driver.fieldAssociation()) << ')' << nl
                << nl
                << "Expression: " << expression << nl
                << "Condition: " << condition << nl
                << nl
                << exit(FatalError);
        }
    }

    if (!ctrl.createNew && driver.resultType() != oldFieldType)
    {
        FatalErrorInFunction
            << "Inconsistent types: " << fieldName << " is  "
            << oldFieldType
            << " but the expression evaluates to "
            << driver.resultType()
            << exit(FatalError);
    }

    Info<< "Dispatch ... " << driver.resultType() << nl;

    #undef setFieldDispatch
    #define setFieldDispatch(FieldType)                                       \
    {                                                                         \
        /* FieldType */                                                       \
        const auto* ptr = driver.isResultType<FieldType>();                   \
        if (ptr)                                                              \
        {                                                                     \
            /* driver.getResult<FieldType>(correctPatches), */                \
                                                                              \
            setField                                                          \
            (                                                                 \
                fieldName,                                                    \
                mesh,                                                         \
                *ptr,                                                         \
                conditionField,                                               \
                dims,                                                         \
                valuePatches,                                                 \
                ctrl                                                          \
            );                                                                \
            return;                                                           \
        }                                                                     \
    }                                                                         \


    setFieldDispatch(volScalarField);
    setFieldDispatch(volVectorField);
    setFieldDispatch(volTensorField);
    setFieldDispatch(volSymmTensorField);
    setFieldDispatch(volSphericalTensorField);

    setFieldDispatch(surfaceScalarField);
    setFieldDispatch(surfaceVectorField);
    setFieldDispatch(surfaceTensorField);
    setFieldDispatch(surfaceSymmTensorField);
    setFieldDispatch(surfaceSphericalTensorField);

    #undef setFieldDispatch
    #define setFieldDispatch(FieldType)                                       \
    {                                                                         \
        /* FieldType */                                                       \
        const auto* ptr = driver.isResultType<FieldType>();                   \
                                                                              \
        if (ptr)                                                              \
        {                                                                     \
            /* driver.getResult<FieldType>(correctPatches), */                \
                                                                              \
            setField                                                          \
            (                                                                 \
                fieldName,                                                    \
                pointMesh::New(mesh),                                         \
                *ptr,                                                         \
                conditionField,                                               \
                dims,                                                         \
                valuePatches,                                                 \
                ctrl                                                          \
            );                                                                \
            return;                                                           \
        }                                                                     \
    }                                                                         \

    setFieldDispatch(pointScalarField);
    setFieldDispatch(pointVectorField);
    setFieldDispatch(pointTensorField);
    setFieldDispatch(pointSymmTensorField);
    setFieldDispatch(pointSphericalTensorField);

    #undef setFieldDispatch

    // Nothing dispatched?

    FatalErrorInFunction
        << "Expression evaluates to an unsupported type: "
        << driver.resultType() << nl << nl
        << "Expression " << expression << nl << endl
        << exit(FatalError);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noFunctionObjects(true);

    // No -constant, no special treatment for 0/
    timeSelector::addOptions(false);

    argList::addOption
    (
        "dict",
        "file",
        "Alternative dictionary for setExprFieldsDict"
    );

    argList::addBoolOption
    (
        "dry-run",
        "Evaluate but do not write"
    );

    argList::addBoolOption
    (
        "verbose",
        "Additional verbosity",
        true // Advanced option
    );

    argList::addOption
    (
        "field",
        "name",
        "The field to overwrite command-line operation)",
        true // Advanced option
    );
    argList::addOption
    (
        "expression",
        "expr",
        "The expression to evaluate (command-line operation)",
        true // Advanced option
    );
    argList::addOption
    (
        "condition",
        "logic",
        "The logical condition when to apply the expression"
        " (command-line operation)",
        true // Advanced option
    );
    argList::addOption
    (
        "dimension",
        "dims",
        "The dimensions to apply for created fields"
        " (command-line operation)",
        true // Advanced option
    );
    argList::addBoolOption
    (
        "debug-parser",
        "Additional debugging information",
        true // Advanced option
    );
    argList::addBoolOption
    (
        "no-variable-cache",
        "Disable caching of expression variables",
        true // Advanced option
    );
    argList::addBoolOption
    (
        "create",
        "Create a new field (command-line operation)",
        true // Advanced option
    );
    argList::addBoolOption
    (
        "keepPatches",
        "Leave patches unaltered"
        " (command-line operation)",
        true // Advanced option
    );
    argList::addOption
    (
        "value-patches",
        "(patches)",
        "A list of patches that receive a fixed value"
        " (command-line operation)",
        true // Advanced option
    );
    argList::addBoolOption
    (
        "dummy-phi",
        "(command-line operation)",
        true // Advanced option
    );

    // Future?
    #if 0
    argList::addBoolOption
    (
        "noCorrectPatches",
        ""
    );
    argList::addBoolOption
    (
        "correctResultBoundaryFields",
        "",
        true
    );
    #endif

    #include "addRegionOption.H"
    #include "setRootCase.H"

    #include "createTime.H"

    const bool dryrun = args.found("dry-run");
    const bool verbose = args.found("verbose");

    const word dictName("setExprFieldsDict");

    instantList times = timeSelector::select0(runTime, args);

    if (times.size() < 1)
    {
        FatalErrorInFunction
            << "No times selected." << exit(FatalError);
    }

    // Disable dimension checking during operation
    dimensionSet::debug = false;

    #include "createNamedMesh.H"

    autoPtr<surfaceScalarField> dummyPhi;

    autoPtr<IOdictionary> exprDictPtr;

    // Sort out conflicts

    const bool useCommandArgs = args.found("field");

    if (useCommandArgs)
    {
        if (args.found("dict"))
        {
            FatalErrorInFunction
                << "Cannot specify both dictionary and command-line arguments"
                << nl
                << endl;
        }

        if (args.found("create") && args.found("keepPatches"))
        {
            FatalErrorInFunction
                << "Cannot specify both 'create' and 'keepPatches'" << nl
                << endl;
        }
    }
    else
    {
        // Carp about inapplicable options

        wordHashSet badOptions
        ({
            "create", "keepPatches", "valuePatches",
            "dimension", "condition", "expression"
        });

        badOptions.retain(args.options());

        if (!badOptions.empty())
        {
            FatalErrorInFunction
                << "Using a dictionary. Cannot specify command options:" << nl
                << nl
                << flatOutput(badOptions.sortedToc()) << nl
                << endl;
        }

        #include "setSystemMeshDictionaryIO.H"
        exprDictPtr.reset(new IOdictionary(dictIO));
    }


    forAll(times, timei)
    {
        runTime.setTime(times[timei], timei);

        Info<< "\nTime = " << runTime.timeName() << endl;

        mesh.readUpdate();

        if (args.found("dummy-phi") && !dummyPhi.valid())
        {
            Info<< "Adding a dummy phi" << endl;
            dummyPhi.reset
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        "phi",
                        mesh.thisDb().time().constant(),
                        mesh.thisDb(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar(Zero)
                )
            );
        }

        if (args.found("withFunctionObjects"))
        {
            runTime.functionObjects().start();
        }

        if (args.found("field"))
        {
            const word fieldName(args.get<word>("field"));

            Info<< "Using command-line options for "
                << fieldName << nl << endl;

            setExprFieldsControl ctrl;

            ctrl.dryRun = dryrun;
            ctrl.debugParsing = args.found("debug-parser");
            ctrl.cacheVariables = !args.found("no-variable-caching");

            ctrl.createNew = args.found("create");
            ctrl.keepPatches = args.found("keepPatches");
            ctrl.correctPatches = !args.found("noCorrectPatches");
            ctrl.correctBCs = args.found("correctResultBoundaryFields");
            ctrl.useDimension = args.found("dimension");

            expressions::exprString
                expression
                (
                    args.opt("expression"),
                    dictionary::null
                );

            expressions::exprString condition;
            if (args.found("condition"))
            {
                args.readIfPresent("condition", condition);
            }

            dimensionSet dims;

            if (ctrl.useDimension)
            {
                ITstream is(args.lookup("dimension"));
                is >> dims;
            }

            evaluate
            (
                mesh,
                fieldName,
                expression,
                condition,
                dictionary::null,
                dims,
                args.getList<word>("valuePatches", false),

                ctrl
            );
        }
        else if (exprDictPtr.valid())
        {
            const dictionary& exprDict = exprDictPtr();

            // Read set construct info from dictionary
            PtrList<entry> actions(exprDict.lookup("expressions"));

            for (const entry& dEntry : actions)
            {
                if (!dEntry.isDict())
                {
                    Info<< "Ignore non-dictionary entry: "
                        << dEntry.keyword() << nl;
                    continue;
                }

                const dictionary& dict = dEntry.dict();

                setExprFieldsControl ctrl;

                ctrl.dryRun = dryrun;
                ctrl.debugParsing = args.found("debug-parser");
                ctrl.cacheVariables = !args.found("no-variable-caching");

                ctrl.createNew = dict.getOrDefault("create", false);
                ctrl.keepPatches = dict.getOrDefault("keepPatches", false);
                ctrl.correctPatches = !args.found("noCorrectPatches");
                ctrl.correctBCs = args.found("correctResultBoundaryFields");

                if (ctrl.createNew && ctrl.keepPatches)
                {
                    FatalIOErrorInFunction(dict)
                        << "Cannot specify both 'create' and 'keepPatches'"
                        << nl << endl
                        << exit(FatalIOError);
                }

                // Local override
                dict.readIfPresent
                (
                    "correctResultBoundaryFields",
                    ctrl.correctBCs
                );


                const word fieldName(dict.get<word>("field"));

                expressions::exprString expression
                (
                    dict.get<string>("expression"),
                    dict
                );

                expressions::exprString condition;

                if (dict.found("condition"))
                {
                    condition =
                        expressions::exprString
                        (
                            dict.get<string>("condition"),
                            dict
                        );
                }

                ctrl.useDimension = dict.found("dimension");

                dimensionSet dims;
                if (ctrl.useDimension)
                {
                    dict.lookup("dimension") >> dims;
                }

                wordList valuePatches;
                dict.readIfPresent("valuePatches", valuePatches);

                if (verbose && !timei)
                {
                    // Report once
                    Info<< "Processing" << dict << nl;
                }

                evaluate
                (
                    mesh,
                    fieldName,
                    expression,
                    condition,
                    dict,
                    dims,
                    valuePatches,

                    ctrl
                );
            }
        }
        else if (exprDictPtr.valid())
        {
            FatalErrorInFunction
                << "No command-line or dictionary??" << nl << endl
                << exit(FatalError);
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
