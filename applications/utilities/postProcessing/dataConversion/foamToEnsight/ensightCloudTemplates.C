/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "ensightCloud.H"
#include "ensightFile.H"
#include "Time.H"
#include "IOField.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "ensightPTraits.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::writeCloudField
(
    const Foam::IOField<Type>& field,
    Foam::ensightFile& os
)
{
    if (returnReduce(field.size(), sumOp<label>()) > 0)
    {
        if (Pstream::master())
        {
            // 6 values per line
            label count = 0;

            // Master
            forAll(field, i)
            {
                Type val = field[i];

                if (mag(val) < 1e-90)
                {
                    val = Zero;
                }

                for (direction i=0; i < pTraits<Type>::nComponents; ++i)
                {
                    label cmpt = ensightPTraits<Type>::componentOrder[i];
                    os.write( component(val, cmpt) );

                    if (++count % 6 == 0)
                    {
                        os.newline();
                    }
                }
            }

            // Slaves
            for (int slave=1; slave<Pstream::nProcs(); ++slave)
            {
                IPstream fromSlave(Pstream::scheduled, slave);
                Field<Type> slaveData(fromSlave);

                forAll(slaveData, i)
                {
                    Type val = slaveData[i];

                    if (mag(val) < 1e-90)
                    {
                        val = Zero;
                    }

                    for (direction i=0; i < pTraits<Type>::nComponents; ++i)
                    {
                        label cmpt = ensightPTraits<Type>::componentOrder[i];
                        os.write( component(val, cmpt) );

                        if (++count % 6 == 0)
                        {
                            os.newline();
                        }
                    }
                }
            }

            // add final newline if required
            if (count % 6)
            {
                os.newline();
            }
        }
        else
        {
            OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
            toMaster<< field;
        }
    }
}


template<class Type>
void Foam::ensightCloudField
(
    const Foam::IOobject& fieldObject,
    const Foam::fileName& dataDir,
    const Foam::label timeIndex,
    const Foam::word& cloudName,
    const Foam::label cloudNo,
    Foam::Ostream& ensightCaseFile,
    const bool dataExists,
    Foam::IOstream::streamFormat format
)
{
    const ensight::VarName varName(fieldObject.name());

    if (dataExists)
    {
        Info<< ' ' << fieldObject.name();
    }
    else
    {
        Info<< ' ' << fieldObject.name() << "{0}"; // ie, empty field
    }

    ensightFile* filePtr(nullptr);
    if (Pstream::master())
    {
        const fileName postFileName =
            ensightFile::subDir(timeIndex)/cloud::prefix/cloudName/varName;

        // the ITER/lagrangian subdirectory must exist
        // the ITER/lagrangian subdirectory was already created
        // when writing positions

        mkDir(dataDir/postFileName.path());

        if (timeIndex == 0)
        {
            const fileName dirName =
                dataDir.name()/ensightFile::mask()/cloud::prefix/cloudName;

            ensightCaseFile.setf(ios_base::left);

            // prefix variables with 'c' (cloud)
            ensightCaseFile
                << ensightPTraits<Type>::typeName << " per "
                << setw(20)
                << "measured node:"
                << " 1  "
                << setw(15)
                << ("c" + Foam::name(cloudNo) + varName).c_str() << ' '
                << (dirName/varName).c_str()
                << nl;
        }

        filePtr = new ensightFile(dataDir, postFileName, format);
        filePtr->write
        (
            // description
            string(postFileName + " <" + pTraits<Type>::typeName + ">")
        );
        filePtr->newline();
    }

    if (dataExists)
    {
        IOField<Type> field(fieldObject);
        writeCloudField(field, *filePtr);
    }

    if (filePtr) // on master only
    {
        delete filePtr;
    }
}


// ************************************************************************* //
