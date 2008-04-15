forAll(Times, n1)
{
    IOobject fieldObjectHeader
    (
        fieldName,
        Times[n1].name(),
        "lagrangian",
        mesh,
        IOobject::NO_READ
    );

    variableGood = variableGood && fieldObjectHeader.headerOk();
}
