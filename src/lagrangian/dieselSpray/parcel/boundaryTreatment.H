if (isA<wallPolyPatch>(pbMesh[patch(face())]))
{
    td.keepParticle = sDB.wall().wallTreatment(*this, face());

    if (sDB.twoD())
    {
        scalar vs = U() & sDB.axisOfSymmetry();
        vector Us = U() - (vs * sDB.axisOfSymmetry());
        scalar magUs = mag(Us);
        scalar v1 = Us & n();
        scalar magV1 = mag(v1);
        U() = (magUs/magV1)*v1*n() + vs*sDB.axisOfSymmetry();
    }
}
else if (isA<wedgePolyPatch>(pbMesh[patch(face())]))
{
    // check if parcel is trying to move out of the domain

    vector nf = normal();

    scalar Un = U() & nf;

    if (Un > 0)
    {
        scalar Un2 = U() & n();
        U() -= 2.0*Un2*n();
    }
}
else if (isA<symmetryPolyPatch>(pbMesh[patch(face())]))
{
    // check if parcel is trying to move out of the domain

    vector nf = normal();

    scalar Un = U() & nf;

    if (Un > 0)
    {
        if (sDB.twoD())
        {
            U() -= 2.0*(U() & n())*n();
        }
        else
        {
            nf /= mag(nf);
            U() -= 2.0*(U() & nf)*nf;
        }
    }
}
else if (isType<polyPatch>(pbMesh[patch(face())]))
{
    // Parcel has hit an inlet or outlet or some such so throw it away
    td.keepParticle = false;
}
