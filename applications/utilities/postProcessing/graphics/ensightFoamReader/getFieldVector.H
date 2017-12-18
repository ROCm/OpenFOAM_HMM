if (nVar >= Num_variables - nSprayVariables)
{
    return Z_UNDEF;
}


IOobject fieldObjectPtr
(
    fieldNames[var2field[nVar]],
    runTime.timeName(),
    mesh,
    IOobject::NO_READ
);

if (!fieldObjectPtr.typeHeaderOk<volVectorField>(true))
{
    return Z_UNDEF;
}

IOobject fieldObject
(
    fieldNames[var2field[nVar]],
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
);

volVectorField vectorField
(
    fieldObject,
    mesh
);

const cellShapeList& cellShapes = meshPtr->cellShapes();

// hexa's
if (which_type == Z_HEX08)
{
    const cellModel& hex = cellModel::ref(cellModel::HEX);
    //const cellModel& wedge = cellModel::ref(cellModel::WEDGE);

    label counter = 1;
    for (label n=0; n<nCells; n++)
    {
        const cellShape& cellShape = cellShapes[n];
        const cellModel& cellModel = cellShape.model();

        if (cellModel == hex) // || (cellModel == wedge))
        {
            var_array[counter++] = vectorField[n][component];
        }
    }
}

// penta's
if (which_type == Z_PEN06)
{
    const cellModel& prism = cellModel::ref(cellModel::PRISM);

    label counter = 1;
    for (label n=0; n<nCells; n++)
    {
        const cellShape& cellShape = cellShapes[n];
        const cellModel& cellModel = cellShape.model();

        if (cellModel == prism)
        {
            var_array[counter++] = vectorField[n][component];
        }
    }
}

// pyramids's
if (which_type == Z_PYR05)
{
    const cellModel& pyr = cellModel::ref(cellModel::PYR);

    label counter = 1;
    for (label n=0; n<nCells; n++)
    {
        const cellShape& cellShape = cellShapes[n];
        const cellModel& cellModel = cellShape.model();

        if (cellModel == pyr)
        {
            var_array[counter++] = vectorField[n][component];
        }
    }
}


// tet's
if (which_type == Z_TET04)
{
    const cellModel& tet = cellModel::ref(cellModel::TET);

    label counter = 1;
    for (label n=0; n<nCells; n++)
    {
        const cellShape& cellShape = cellShapes[n];
        const cellModel& cellModel = cellShape.model();

        if (cellModel == tet)
        {
            var_array[counter++] = vectorField[n][component];
        }
    }
}

if (which_type == Z_NFACED)
{
    const cellList& cells = meshPtr->cells();

    label counter = 1;
    for (label n=0; n<nCells; n++)
    {
        const labelList& points = cellShapes[n];
        label nFacesInCell = cells[n].size();

        if ((nFacesInCell == 6) && (points.size() == 8))
        {}
        else if ((nFacesInCell == 4) && (points.size() == 4))
        {}
        else if (nFacesInCell == 5)
        {
            if (points.size() == 6)
            {}
            else if (points.size() == 5)
            {}
            else
            {
                var_array[counter++] = vectorField[n][component];
            }
        }
        else
        {
            var_array[counter++] = vectorField[n][component];
        }
    }
}
