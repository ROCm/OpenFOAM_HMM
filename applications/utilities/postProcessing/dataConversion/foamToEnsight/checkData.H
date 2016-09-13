// ignore special fields or fields that we don't handle
//
bool variableGood = true;
for (label n1=0; n1<timeDirs.size() && variableGood; ++n1)
{
    // ignore _0 fields
    if (fieldName.size() > 2 && fieldName(fieldName.size() - 2, 2) == "_0")
    {
        variableGood = false;
    }
    else
    {
        variableGood = IOobject
        (
            fieldName,
            timeDirs[n1].name(),
            mesh,
            IOobject::NO_READ
        ).typeHeaderOk<volScalarField>(false);
    }
}
