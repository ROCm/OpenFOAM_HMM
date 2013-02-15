/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "lduPrimitiveMesh.H"
#include "processorLduInterface.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lduPrimitiveMesh::lduPrimitiveMesh
(
    const label nCells,
    const labelUList& l,
    const labelUList& u,
    const labelListList& pa,
    const lduInterfacePtrsList& interfaces,
    const lduSchedule& ps,
    const label comm
)
:
    lduAddressing(nCells),
    lowerAddr_(l),
    upperAddr_(u),
    patchAddr_(pa),
    interfaces_(interfaces),
    patchSchedule_(ps),
    comm_(comm)
{
    Pout<< "lduPrimitiveMesh :"
        << " nCells:" << nCells
        << " l:" << lowerAddr_.size()
        << " u:" << upperAddr_.size()
        << " pa:" << patchAddr_.size()
        << " interfaces:" << interfaces_.size()
        << " comm:" << comm_
        << endl;
    forAll(interfaces_, i)
    {
        if (interfaces_.set(i))
        {
            if (isA<processorLduInterface>(interfaces_[i]))
            {
                const processorLduInterface& pi = refCast
                <
                    const processorLduInterface
                >(interfaces_[i]);

                Pout<< "    patch:" << i
                    << " size:" << patchAddr_[i].size()
                    << " myProcNo:" << pi.myProcNo()
                    << " neighbProcNo:" << pi.neighbProcNo()
                    << " comm:" << pi.comm()
                    << endl;
            }
        }
    }
}


Foam::lduPrimitiveMesh::lduPrimitiveMesh
(
    const label nCells,
    labelList& l,
    labelList& u,
    labelListList& pa,
    lduInterfacePtrsList& interfaces,
    const lduSchedule& ps,
    const label comm,
    bool reUse
)
:
    lduAddressing(nCells),
    lowerAddr_(l, reUse),
    upperAddr_(u, reUse),
    patchAddr_(pa, reUse),
    interfaces_(interfaces, reUse),
    patchSchedule_(ps),
    comm_(comm)
{
    Pout<< "lduPrimitiveMesh :"
        << " nCells:" << nCells
        << " l:" << lowerAddr_.size()
        << " u:" << upperAddr_.size()
        << " pa:" << patchAddr_.size()
        << " interfaces:" << interfaces_.size()
        << " comm:" << comm_
        << endl;
    forAll(interfaces_, i)
    {
        if (interfaces_.set(i))
        {
            if (isA<processorLduInterface>(interfaces_[i]))
            {
                const processorLduInterface& pi = refCast
                <
                    const processorLduInterface
                >(interfaces_[i]);

                Pout<< "    patch:" << i
                    << " size:" << patchAddr_[i].size()
                    << " myProcNo:" << pi.myProcNo()
                    << " neighbProcNo:" << pi.neighbProcNo()
                    << " comm:" << pi.comm()
                    << endl;
            }
        }
    }
}


//Foam::lduPrimitiveMesh::lduPrimitiveMesh(Istream& is)
//:
//    lduAddressing(readLabel(is)),
//    lowerAddr_(is),
//    upperAddr_(is),
//    patchAddr_(is),
//    interfaces_(is),
//    patchSchedule_(is),
//    comm_(readLabel(is))
//{}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//Foam::Ostream& Foam::operator<<
//(
//    Ostream& os,
//    const Foam::lduPrimitiveMesh& mesh
//)
//{
//    os  << mesh.size()
//        << mesh.lowerAddr_
//        << mesh.upperAddr_
//        << mesh.patchAddr_
//        << mesh.interfaces_
//        << mesh.patchSchedule_
//        << mesh.comm_
//
//    return os;
//}


// ************************************************************************* //
