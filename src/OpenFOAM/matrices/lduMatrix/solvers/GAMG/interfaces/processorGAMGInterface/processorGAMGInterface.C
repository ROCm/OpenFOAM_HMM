/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "processorGAMGInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "Map.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorGAMGInterface, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterface,
        processorGAMGInterface,
        lduInterface
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorGAMGInterface::processorGAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const labelField& localRestrictAddressing,
    const labelField& neighbourRestrictAddressing
)
:
    GAMGInterface
    (
        index,
        coarseInterfaces,
        fineInterface,
        localRestrictAddressing,
        neighbourRestrictAddressing
    ),
    fineProcInterface_(refCast<const processorLduInterface>(fineInterface))
{
    // Make a lookup table of entries for owner/neighbour
    Map<SLList<label> > neighboursTable
    (
        localRestrictAddressing.size()
    );

    // Table of face-sets to be agglomerated
    Map<SLList<SLList<label> > > faceFaceTable
    (
        localRestrictAddressing.size()
    );

    label nCoarseFaces = 0;

    forAll(localRestrictAddressing, ffi)
    {
        label curMaster = -1;
        label curSlave = -1;

        // Do switching on master/slave indexes based on the owner/neighbour of
        // the processor index such that both sides get the same answer.
        if (myProcNo() < neighbProcNo())
        {
            // Master side
            curMaster = localRestrictAddressing[ffi];
            curSlave = neighbourRestrictAddressing[ffi];
        }
        else
        {
            // Slave side
            curMaster = neighbourRestrictAddressing[ffi];
            curSlave = localRestrictAddressing[ffi];
        }

        // Look for the master cell.  If it has already got a face,
        // add the coefficient to the face.  If not, create a new face.
        if (neighboursTable.found(curMaster))
        {
            // Check all current neighbours to see if the current slave already
            // exists and if so, add the fine face to the agglomeration.

            SLList<label>& curNbrs = neighboursTable.find(curMaster)();

            SLList<SLList<label> >& curFaceFaces =
                faceFaceTable.find(curMaster)();

            bool nbrFound = false;

            SLList<label>::iterator nbrsIter = curNbrs.begin();

            SLList<SLList<label> >::iterator faceFacesIter =
                curFaceFaces.begin();

            for
            (
                ;
                nbrsIter != curNbrs.end(), faceFacesIter != curFaceFaces.end();
                ++nbrsIter, ++faceFacesIter
            )
            {
                if (nbrsIter() == curSlave)
                {
                    nbrFound = true;
                    faceFacesIter().append(ffi);
                    break;
                }
            }

            if (!nbrFound)
            {
                curNbrs.append(curSlave);
                curFaceFaces.append(ffi);

                // New coarse face created
                nCoarseFaces++;
            }
        }
        else
        {
            // This master has got no neighbours yet.  Add a neighbour
            // and a coefficient, thus creating a new face
            neighboursTable.insert(curMaster, SLList<label>(curSlave));
            faceFaceTable.insert(curMaster, SLList<SLList<label> >(ffi));

            // New coarse face created
            nCoarseFaces++;
        }
    } // end for all fine faces



    faceCells_.setSize(nCoarseFaces, -1);
    faceRestrictAddressing_.setSize(localRestrictAddressing.size());

    labelList contents = neighboursTable.toc();

    // Reset face counter for re-use
    nCoarseFaces = 0;

    if (myProcNo() < neighbProcNo())
    {
        // On master side, the owner addressing is stored in table of contents
        forAll(contents, masterI)
        {
            SLList<label>& curNbrs = neighboursTable.find(contents[masterI])();

            SLList<SLList<label> >& curFaceFaces =
                faceFaceTable.find(contents[masterI])();

            SLList<label>::iterator nbrsIter = curNbrs.begin();

            SLList<SLList<label> >::iterator faceFacesIter =
                curFaceFaces.begin();

            for
            (
                ;
                nbrsIter != curNbrs.end(), faceFacesIter != curFaceFaces.end();
                ++nbrsIter, ++faceFacesIter
            )
            {
                faceCells_[nCoarseFaces] = contents[masterI];

                forAllConstIter(SLList<label>, faceFacesIter(), facesIter)
                {
                    faceRestrictAddressing_[facesIter()] = nCoarseFaces;
                }

                nCoarseFaces++;
            }
        }
    }
    else
    {
        // On slave side, the owner addressing is stored in linked lists
        forAll(contents, masterI)
        {
            SLList<label>& curNbrs = neighboursTable.find(contents[masterI])();

            SLList<SLList<label> >& curFaceFaces =
                faceFaceTable.find(contents[masterI])();

            SLList<label>::iterator nbrsIter = curNbrs.begin();

            SLList<SLList<label> >::iterator faceFacesIter =
                curFaceFaces.begin();

            for
            (
                ;
                nbrsIter != curNbrs.end(), faceFacesIter != curFaceFaces.end();
                ++nbrsIter, ++faceFacesIter
            )
            {
                faceCells_[nCoarseFaces] = nbrsIter();

                forAllConstIter(SLList<label>, faceFacesIter(), facesIter)
                {
                    faceRestrictAddressing_[facesIter()] = nCoarseFaces;
                }

                nCoarseFaces++;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Desstructor * * * * * * * * * * * * * * * //

Foam::processorGAMGInterface::~processorGAMGInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::processorGAMGInterface::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    send(commsType, interfaceInternalField(iF)());
}


Foam::tmp<Foam::labelField> Foam::processorGAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    return receive<label>(commsType, this->size());
}


// ************************************************************************* //
