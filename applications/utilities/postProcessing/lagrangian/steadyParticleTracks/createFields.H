const word dictName("particleTrackDict");

#include "setConstantMeshDictionaryIO.H"

IOdictionary propsDict(dictIO);

word cloudName(propsDict.lookup("cloud"));

List<word> userFields(propsDict.lookup("fields"));
