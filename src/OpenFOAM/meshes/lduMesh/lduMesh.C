/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "lduMesh.H"
#include "objectRegistry.H"
#include "processorLduInterface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(lduMesh, 0);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::objectRegistry& Foam::lduMesh::thisDb() const
{
    notImplemented("lduMesh::thisDb() const");
    const objectRegistry* orPtr_ = NULL;
    return *orPtr_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const InfoProxy<lduMesh>& ip)
{
    const lduMesh& ldum = ip.t_;
    const lduAddressing& addr = ldum.lduAddr();
    const lduInterfacePtrsList interfaces = ldum.interfaces();

    Pout<< "lduMesh :"
        << " size:" << addr.size()
        << " l:" << addr.lowerAddr().size()
        << " u:" << addr.upperAddr().size()
        << " interfaces:" << interfaces.size()
        << " comm:" << ldum.comm()
        << endl;
    forAll(interfaces, i)
    {
        if (interfaces.set(i))
        {
            const labelUList& faceCells = addr.patchAddr(i);

            if (isA<processorLduInterface>(interfaces[i]))
            {
                const processorLduInterface& pi = refCast
                <
                    const processorLduInterface
                >(interfaces[i]);

                Pout<< "    patch:" << i
                    << " type:" << interfaces[i].type()
                    << " size:" << faceCells.size()
                    << " myProcNo:" << pi.myProcNo()
                    << " neighbProcNo:" << pi.neighbProcNo()
                    << " comm:" << pi.comm()
                    << endl;
            }
            else
            {
                Pout<< "    patch:" << i
                    << " type:" << interfaces[i].type()
                    << " size:" << faceCells.size()
                    << endl;
            }
        }
    }


    // Print actual contents
    {
        const labelList& l = addr.lowerAddr();
        const labelList& u = addr.upperAddr();
        forAll(l, faceI)
        {
            Pout<< "        face:" << faceI << " l:" << l[faceI]
                << " u:" << u[faceI] << endl;
        }
        forAll(interfaces, i)
        {
            if (interfaces.set(i))
            {
                const labelUList& faceCells = addr.patchAddr(i);
                if (faceCells.size())
                {
                    Pout<< "    patch:" << i
                        << " type:" << interfaces[i].type() << endl;

                    if (isA<processorLduInterface>(interfaces[i]))
                    {
                        const processorLduInterface& pi = refCast
                        <
                            const processorLduInterface
                        >(interfaces[i]);

                        Pout<< "    myProcNo:" << pi.myProcNo()
                            << " neighbProcNo:" << pi.neighbProcNo()
                            << " comm:" << pi.comm()
                            << endl;
                    }

                    forAll(faceCells, i)
                    {
                        Pout<< "        " << i << " own:" << faceCells[i]
                            << endl;
                    }
                }
            }
        }
    }

    os.check("Ostream& operator<<(Ostream&, const lduMesh&");

    return os;
}


// ************************************************************************* //
