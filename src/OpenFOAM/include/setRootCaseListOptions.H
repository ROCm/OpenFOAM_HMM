// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Declare some "standard" list options

#ifndef setRootCaseListOptions_H
#define setRootCaseListOptions_H

argList::addBoolOption
(
    "listSwitches",
    "List switches declared in libraries but not set in etc/controlDict",
    true // advanced
);
argList::addBoolOption
(
    "listRegisteredSwitches",
    "List switches registered for run-time modification",
    true // advanced
);
argList::addBoolOption
(
    "listUnsetSwitches",
    "List switches declared in libraries but not set in etc/controlDict",
    true // advanced
);

#ifdef fvPatchField_H
argList::addBoolOption
(
    "listScalarBCs",
    "List scalar field boundary conditions (fvPatchField<scalar>)",
    true // advanced
);
argList::addBoolOption
(
    "listVectorBCs",
    "List vector field boundary conditions (fvPatchField<vector>)",
    true // advanced
);
#endif

#ifdef functionObject_H
argList::addBoolOption
(
    "listFunctionObjects",
    "List functionObjects",
    true // advanced
);
#endif

#ifdef fvOption_H
argList::addBoolOption
(
    "listFvOptions",
    "List fvOptions",
    true // advanced
);
#endif

#if defined(turbulentTransportModel_H) || defined(turbulentFluidThermoModel_H)
argList::addBoolOption
(
    "listTurbulenceModels",
    "List turbulenceModels",
    true // advanced
);
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
