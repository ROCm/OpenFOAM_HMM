// ignore special fields (_0 fields),
// ignore fields we don't handle,
// ignore fields that are not available for all time-steps

// hash by field-name in fieldsToUse
if (!fieldsToUse.found(fieldName))
{
    bool variableGood = false;

    forAll(timeDirs, n1)
    {
        variableGood =
        (
            !fieldName.endsWith("_0")
         && IOobject
            (
                fieldName,
                timeDirs[n1].name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false  // no register
            ).typeHeaderOk<volScalarField>(false, false)
        );

        if (!variableGood)
        {
            break;
        }
    }

    reduce(variableGood, andOp<bool>());
    fieldsToUse.set(fieldName, variableGood);
}
