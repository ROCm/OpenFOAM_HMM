// Read particleTrackDict dictionary and extract values from it

const word dictName("particleTrackDict");

#include "setConstantMeshDictionaryIO.H"

IOdictionary propsDict(dictIO);

const word cloudName(propsDict.get<word>("cloud"));

List<word> userFields(propsDict.lookup("fields"));

const dictionary formatOptions
(
    propsDict.subOrEmptyDict("formatOptions", keyType::LITERAL)
);


// ************************************************************************* //
