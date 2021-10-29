/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Namespace
    Foam::ensightOutput

Description
    A collection of functions for writing ensight file content.

SourceFiles
    ensightOutput.C
    ensightOutputTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef ensightOutput_H
#define ensightOutput_H

#include "ensightFile.H"
#include "ensightGeoFile.H"
#include "ensightCells.H"
#include "ensightFaces.H"
#include "ensightPTraits.H"

#include "faceListFwd.H"
#include "cellListFwd.H"

#include "ListOps.H"
#include "ListListOps.H"
#include "IndirectList.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Local definitions, could be relocated to ListListOps directly

namespace Foam
{
namespace ListListOps
{

//- Return the sizes of the sub-lists
template<class T, class Addr, class AccessOp>
labelList subSizes
(
    const IndirectListBase<T, Addr>& lists,
    AccessOp aop
)
{
    labelList output(lists.size());
    auto out = output.begin();

    for (const T& sub : lists)
    {
        *out = aop(sub).size();
        ++out;
    }

    return output;
}


//- Inplace renumber the values (not the indices) of a list of lists.
//  Negative IntListType elements are left untouched.
template<class IntListType>
void inplaceRenumber
(
    const labelUList& oldToNew,
    IntListType& lists
)
{
    for (auto& sub : lists)
    {
        for (auto& item : sub)
        {
            if (item >= 0)
            {
                item = oldToNew[item];
            }
        }
    }
}

} // End namespace ListListOps
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class cellShape;
class polyMesh;

namespace ensightOutput
{

/*---------------------------------------------------------------------------*\
                                Geometry Output
\*---------------------------------------------------------------------------*/

//- Write list of faces
void writeFaceList
(
    ensightGeoFile& os,
    const UList<face>& faces
);

//- Write list of faces
void writeFaceList
(
    ensightGeoFile& os,
    const UIndirectList<face>& faces
);

//- Write cell connectivity via cell shapes
void writeCellShapes
(
    ensightGeoFile& os,
    const UList<cellShape>& shapes
);


//- Write the point ids per poly element.
//  Points have been already renumbered
void writePolysPoints
(
    ensightGeoFile& os,
    const cellUList& meshCells,
    const labelUList& addr,     //!< Cell ids to write
    const faceUList& meshFaces,
    const labelUList& faceOwner
);

//- Write the point ids per poly element, with point renumbering
void writePolysPoints
(
    ensightGeoFile& os,
    const polyMesh& mesh,
    const labelUList& addr,     //!< Cell ids to write
    const labelList& pointMap   //!< Point map to use
);


//- Write the regular face connectivity for specified type and
//- and specified faces
void writeFaceConnectivity
(
    ensightGeoFile& os,
    const ensightFaces::elemType etype,
    const label nTotal,
    const UIndirectList<face>& faces,
    bool parallel               //!< Collective write?
);


//- Write the regular face connectivity for specified type
void writeFaceConnectivity
(
    ensightGeoFile& os,
    const ensightFaces::elemType etype,
    const label nTotal,
    const faceUList& faces,
    bool parallel               //!< Collective write?
);


//- Write the face connectivity for the part
void writeFaceConnectivity
(
    ensightGeoFile& os,
    const ensightFaces& part,
    const faceUList& faces,
    bool parallel               //!< Collective write?
);


//- Write the \b presorted face connectivity for the part
//  This is a special case when the list of faces is already in
//  ensight sorted order
void writeFaceConnectivityPresorted
(
    ensightGeoFile& os,
    const ensightFaces& part,
    const faceUList& faces,
    bool parallel               //!< Collective write?
);


/*---------------------------------------------------------------------------*\
                                Field Output
\*---------------------------------------------------------------------------*/

//- Write a field of cell values as an indirect list,
//- using the cell ids from ensightCells
template<class Type>
bool writeField
(
    ensightFile& os,
    const Field<Type>& fld,
    const ensightCells& part,
    bool parallel  //!< Collective write?
);

//- Write a field of faces values as an indirect list,
//- using the face ids from ensightFaces
template<class Type>
bool writeField
(
    ensightFile& os,
    const Field<Type>& fld,
    const ensightFaces& part,
    bool parallel  //!< Collective write?
);


/*---------------------------------------------------------------------------*\
                        Namespace ensightOutput::Detail
\*---------------------------------------------------------------------------*/

//- \brief Implementation details and output backends that would not normally
//- be called directly by a user.

namespace Detail
{

//- Return sizes of faces in the list
labelList getFaceSizes(const UList<face>& faces);

//- Return sizes of faces in the list
labelList getFaceSizes(const UIndirectList<face>& faces);

//- The number of faces per poly element
labelList getPolysNFaces(const polyMesh& mesh, const labelUList& addr);

//- The number of points for each face of the poly elements
labelList getPolysNPointsPerFace(const polyMesh& mesh, const labelUList& addr);


//- Copy specified field component into a scalar buffer
//- works for various lists types. Must be adequately sized before calling
template<template<typename> class FieldContainer, class Type>
void copyComponent
(
    List<scalar>& cmptBuffer,
    const FieldContainer<Type>& input,
    const direction cmpt
);


//- Write field content (component-wise)
template<template<typename> class FieldContainer, class Type>
void writeFieldContent
(
    ensightFile& os,
    const FieldContainer<Type>& fld,
    bool parallel  //!< Collective write?
);


//- Write coordinates (component-wise) for the given part
template<template<typename> class FieldContainer>
bool writeCoordinates
(
    ensightGeoFile& os,
    const label partId,
    const word& partName,
    const label nPoints,
    const FieldContainer<Foam::point>& fld,
    bool parallel  //!< Collective write?
);


//- Write field content (component-wise) for the given ensight element type
template<template<typename> class FieldContainer, class Type>
bool writeFieldComponents
(
    ensightFile& os,
    const char* key,
    const FieldContainer<Type>& fld,
    bool parallel  //!< Collective write?
);


//- Write a sub-field of faces values as an indirect list,
//- using the sub-list sizing information from ensightFaces
template<class Type>
bool writeFaceSubField
(
    ensightFile& os,
    const Field<Type>& fld,
    const ensightFaces& part,
    bool parallel  //!< Collective write?
);


//- Write a field of faces values as an indirect list,
//- using the face order from ensightFaces
template<class Type>
bool writeFaceLocalField
(
    ensightFile& os,
    const Field<Type>& fld,
    const ensightFaces& part,
    bool parallel  //!< Collective write?
);


} // End namespace Detail


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ensightOutput
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ensightOutputTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
