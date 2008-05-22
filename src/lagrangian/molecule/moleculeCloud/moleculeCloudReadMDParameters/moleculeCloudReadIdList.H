IOdictionary idListDict
(
    IOobject
    (
        "idList",
        mesh_.time().constant(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

idList_ = List<word>(idListDict.lookup("idList"));

