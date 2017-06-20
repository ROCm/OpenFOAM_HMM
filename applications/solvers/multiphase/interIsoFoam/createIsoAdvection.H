// Defining isoAdvection
isoAdvection advector(alpha1, phi, U);

bool frozenFlow = pimple.dict().lookupOrDefault<bool>("frozenFlow", false);
if (frozenFlow)
{
    Info<< "Employing frozen-flow assumption: "
        << "pressure-velocity system will not be updated"
        << endl;
}
