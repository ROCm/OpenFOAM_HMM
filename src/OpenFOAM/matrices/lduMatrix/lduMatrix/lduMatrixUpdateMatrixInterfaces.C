/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "lduMatrix.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::lduMatrix::initMatrixInterfaces
(
    const bool add,
    const FieldField<Field, scalar>& coupleCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const solveScalarField& psiif,
    solveScalarField& result,
    const direction cmpt
) const
{
    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(interfaces, interfacei)
        {
            if (interfaces.set(interfacei))
            {
                interfaces[interfacei].initInterfaceMatrixUpdate
                (
                    result,
                    add,
                    mesh().lduAddr(),
                    interfacei,
                    psiif,
                    coupleCoeffs[interfacei],
                    cmpt,
                    Pstream::defaultCommsType
                );
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = this->patchSchedule();

        // Loop over the "global" patches are on the list of interfaces but
        // beyond the end of the schedule which only handles "normal" patches
        for
        (
            label interfacei=patchSchedule.size()/2;
            interfacei<interfaces.size();
            interfacei++
        )
        {
            if (interfaces.set(interfacei))
            {
                interfaces[interfacei].initInterfaceMatrixUpdate
                (
                    result,
                    add,
                    mesh().lduAddr(),
                    interfacei,
                    psiif,
                    coupleCoeffs[interfacei],
                    cmpt,
                    Pstream::commsTypes::blocking
                );
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType]
            << exit(FatalError);
    }
}


void Foam::lduMatrix::updateMatrixInterfaces
(
    const bool add,
    const FieldField<Field, scalar>& coupleCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const solveScalarField& psiif,
    solveScalarField& result,
    const direction cmpt,
    const label startRequest
) const
{
    if (Pstream::defaultCommsType == Pstream::commsTypes::blocking)
    {
        forAll(interfaces, interfacei)
        {
            if (interfaces.set(interfacei))
            {
                interfaces[interfacei].updateInterfaceMatrix
                (
                    result,
                    add,
                    mesh().lduAddr(),
                    interfacei,
                    psiif,
                    coupleCoeffs[interfacei],
                    cmpt,
                    Pstream::defaultCommsType
                );
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking)
    {
        // Try and consume interfaces as they become available
        bool allUpdated = false;

        for (label i=0; i<UPstream::nPollProcInterfaces; i++)
        {
            allUpdated = true;

            forAll(interfaces, interfacei)
            {
                if (interfaces.set(interfacei))
                {
                    if (!interfaces[interfacei].updatedMatrix())
                    {
                        if (interfaces[interfacei].ready())
                        {
                            interfaces[interfacei].updateInterfaceMatrix
                            (
                                result,
                                add,
                                mesh().lduAddr(),
                                interfacei,
                                psiif,
                                coupleCoeffs[interfacei],
                                cmpt,
                                Pstream::defaultCommsType
                            );
                        }
                        else
                        {
                            allUpdated = false;
                        }
                    }
                }
            }

            if (allUpdated)
            {
                break;
            }
        }

        // Block for everything
        if (Pstream::parRun())
        {
            if (allUpdated)
            {
                // All received. Just remove all outstanding requests
                UPstream::resetRequests(startRequest);
            }
            else
            {
                // Block for all requests and remove storage
                UPstream::waitRequests(startRequest);
            }
        }

        // Consume
        forAll(interfaces, interfacei)
        {
            if
            (
                interfaces.set(interfacei)
            && !interfaces[interfacei].updatedMatrix()
            )
            {
                interfaces[interfacei].updateInterfaceMatrix
                (
                    result,
                    add,
                    mesh().lduAddr(),
                    interfacei,
                    psiif,
                    coupleCoeffs[interfacei],
                    cmpt,
                    Pstream::defaultCommsType
                );
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = this->patchSchedule();

        // Loop over all the "normal" interfaces relating to standard patches
        for (const auto& sched : patchSchedule)
        {
            const label interfacei = sched.patch;

            if (interfaces.set(interfacei))
            {
                if (sched.init)
                {
                    interfaces[interfacei].initInterfaceMatrixUpdate
                    (
                        result,
                        add,
                        mesh().lduAddr(),
                        interfacei,
                        psiif,
                        coupleCoeffs[interfacei],
                        cmpt,
                        Pstream::commsTypes::scheduled
                    );
                }
                else
                {
                    interfaces[interfacei].updateInterfaceMatrix
                    (
                        result,
                        add,
                        mesh().lduAddr(),
                        interfacei,
                        psiif,
                        coupleCoeffs[interfacei],
                        cmpt,
                        Pstream::commsTypes::scheduled
                    );
                }
            }
        }

        // Loop over the "global" patches are on the list of interfaces but
        // beyond the end of the schedule which only handles "normal" patches
        for
        (
            label interfacei=patchSchedule.size()/2;
            interfacei<interfaces.size();
            interfacei++
        )
        {
            if (interfaces.set(interfacei))
            {
                interfaces[interfacei].updateInterfaceMatrix
                (
                    result,
                    add,
                    mesh().lduAddr(),
                    interfacei,
                    psiif,
                    coupleCoeffs[interfacei],
                    cmpt,
                    Pstream::commsTypes::blocking
                );
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType]
            << exit(FatalError);
    }
}


// ************************************************************************* //
