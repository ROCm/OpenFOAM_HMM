/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify i
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

#include "externalCoupledFunctionObject.H"
#include "OSspecific.H"
#include "IFstream.H"
#include "OFstream.H"
#include "volFields.H"
#include "externalCoupledMixedFvPatchFields.H"
#include "mixedFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "OStringStream.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::externalCoupledFunctionObject::readData
(
    const fvMesh& mesh,
    const wordRe& groupName,
    const labelList& patchIDs,
    const word& fieldName
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;
    typedef externalCoupledMixedFvPatchField<Type> patchFieldType;

    if (!mesh.foundObject<volFieldType>(fieldName))
    {
        return false;
    }

    const volFieldType& cvf = mesh.lookupObject<volFieldType>(fieldName);
    const typename volFieldType::GeometricBoundaryField& bf =
        cvf.boundaryField();


    // File only opened on master; contains data for all processors, for all
    // patchIDs.
    autoPtr<IFstream> masterFilePtr;
    if (Pstream::master())
    {
        const fileName transferFile
        (
            groupDir(commsDir_, mesh.dbDir(), groupName)
          / fieldName + ".in"
        );

        if (log_)
        {
            Info<< type() << ": reading data from " << transferFile << endl;
        }
        masterFilePtr.reset(new IFstream(transferFile));

        if (!masterFilePtr().good())
        {
            FatalIOErrorIn
            (
                "void externalCoupledFunctionObject::readData"
                "("
                    "const fvMesh&, "
                    "const wordRe&, "
                    "const labelList&, "
                    "const word&"
               ")",
                masterFilePtr()
            )   << "Cannot open file for region " << mesh.name()
                << ", field " << fieldName << ", patches " << patchIDs
                << exit(FatalIOError);
        }
    }

    // Handle column-wise reading of patch data. Supports most easy types
    forAll(patchIDs, i)
    {
        label patchI = patchIDs[i];

        if (isA<patchFieldType>(bf[patchI]))
        {
            // Explicit handling of externalCoupledObjectMixed bcs - they have
            // specialised reading routines.

            patchFieldType& pf = const_cast<patchFieldType&>
            (
                refCast<const patchFieldType>
                (
                    bf[patchI]
                )
            );

            // Read from master into local stream
            OStringStream os;
            readLines
            (
                bf[patchI].size(),      // number of lines to read
                masterFilePtr,
                os
            );

            // Pass responsability for all reading over to bc
            pf.readData(IStringStream(os.str())());

            // Update the value from the read coefficicient. Bypass any
            // additional processing by derived type.
            pf.patchFieldType::evaluate();
        }
        else if (isA<mixedFvPatchField<Type> >(bf[patchI]))
        {
            // Read columns from file for
            // value, snGrad, refValue, refGrad, valueFraction
            List<scalarField> data;
            readColumns
            (
                bf[patchI].size(),              // number of lines to read
                4*pTraits<Type>::nComponents+1, // nColumns: 4*Type + 1*scalar
                masterFilePtr,
                data
            );

            mixedFvPatchField<Type>& pf = const_cast<mixedFvPatchField<Type>&>
            (
                refCast<const mixedFvPatchField<Type> >
                (
                    bf[patchI]
                )
            );

            // Transfer read data to bc.
            // Skip value, snGrad
            direction columnI = 2*pTraits<Type>::nComponents;

            Field<Type>& refValue = pf.refValue();
            for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
            {
                refValue.replace(cmpt, data[columnI++]);
            }
            Field<Type>& refGrad = pf.refGrad();
            for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
            {
                refGrad.replace(cmpt, data[columnI++]);
            }
            pf.valueFraction() = data[columnI];

            // Update the value from the read coefficicient. Bypass any
            // additional processing by derived type.
            pf.mixedFvPatchField<Type>::evaluate();
        }
        else if (isA<fixedGradientFvPatchField<Type> >(bf[patchI]))
        {
            // Read columns for value and gradient
            List<scalarField> data;
            readColumns
            (
                bf[patchI].size(),              // number of lines to read
                2*pTraits<Type>::nComponents,   // nColumns: Type
                masterFilePtr,
                data
            );

            fixedGradientFvPatchField<Type>& pf =
            const_cast<fixedGradientFvPatchField<Type>&>
            (
                refCast<const fixedGradientFvPatchField<Type> >
                (
                    bf[patchI]
                )
            );

            // Transfer gradient to bc
            Field<Type>& gradient = pf.gradient();
            for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
            {
                gradient.replace(cmpt, data[pTraits<Type>::nComponents+cmpt]);
            }

            // Update the value from the read coefficicient. Bypass any
            // additional processing by derived type.
            pf.fixedGradientFvPatchField<Type>::evaluate();
        }
        else if (isA<fixedValueFvPatchField<Type> >(bf[patchI]))
        {
            // Read columns for value only
            List<scalarField> data;
            readColumns
            (
                bf[patchI].size(),              // number of lines to read
                pTraits<Type>::nComponents,     // number of columns to read
                masterFilePtr,
                data
            );

            // Transfer read value to bc
            Field<Type> value(bf[patchI].size());
            for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
            {
                value.replace(cmpt, data[cmpt]);
            }

            fixedValueFvPatchField<Type>& pf =
            const_cast<fixedValueFvPatchField<Type>&>
            (
                refCast<const fixedValueFvPatchField<Type> >
                (
                    bf[patchI]
                )
            );

            pf == value;

            // Update the value from the read coefficicient. Bypass any
            // additional processing by derived type.
            pf.fixedValueFvPatchField<Type>::evaluate();
        }
        else
        {
            FatalErrorIn
            (
                "void externalCoupledFunctionObject::readData"
                "("
                    "const fvMesh&, "
                    "const wordRe&, "
                    "const labelList&, "
                    "const word&"
               ")"
            )
                << "Unsupported boundary condition " << bf[patchI].type()
                << " for patch " << bf[patchI].patch().name()
                << " in region " << mesh.name()
                << exit(FatalError);
        }

        initialised_ = true;
    }

    return true;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::externalCoupledFunctionObject::gatherAndCombine
(
    const Field<Type>& fld
)
{
    // Collect values from all processors
    List<Field<Type> > gatheredValues(Pstream::nProcs());
    gatheredValues[Pstream::myProcNo()] = fld;
    Pstream::gatherList(gatheredValues);


    tmp<Field<Type> > tresult(new Field<Type>(0));
    Field<Type>& result = tresult();

    if (Pstream::master())
    {
        // Combine values into single field
        label globalElemI = 0;

        forAll(gatheredValues, lstI)
        {
            globalElemI += gatheredValues[lstI].size();
        }

        result.setSize(globalElemI);

        globalElemI = 0;

        forAll(gatheredValues, lstI)
        {
            const Field<Type>& sub = gatheredValues[lstI];

            forAll(sub, elemI)
            {
                result[globalElemI++] = sub[elemI];
            }
        }
    }

    return tresult;
}


template<class Type>
bool Foam::externalCoupledFunctionObject::writeData
(
    const fvMesh& mesh,
    const wordRe& groupName,
    const labelList& patchIDs,
    const word& fieldName
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;
    typedef externalCoupledMixedFvPatchField<Type> patchFieldType;

    if (!mesh.foundObject<volFieldType>(fieldName))
    {
        return false;
    }

    const volFieldType& cvf = mesh.lookupObject<volFieldType>(fieldName);
    const typename volFieldType::GeometricBoundaryField& bf =
        cvf.boundaryField();


    // File only opened on master; contains data for all processors, for all
    // patchIDs
    autoPtr<OFstream> masterFilePtr;
    if (Pstream::master())
    {
        const fileName transferFile
        (
            groupDir(commsDir_, mesh.dbDir(), groupName)
          / fieldName + ".out"
        );

        if (log_)
        {
            Info<< type() << ": writing data to " << transferFile << endl;
        }
        masterFilePtr.reset(new OFstream(transferFile));

        if (!masterFilePtr().good())
        {
            FatalIOErrorIn
            (
                "externalCoupledFunctionObject::writeData"
                "("
                    "const fvMesh&, "
                    "const wordRe&, "
                    "const labelList&, "
                    "const word&"
                ") const",
                masterFilePtr()
            )   << "Cannot open file for region " << mesh.name()
                << ", field " << fieldName << ", patches " << patchIDs
                << exit(FatalIOError);
        }
    }


    bool headerDone = false;

    // Handle column-wise writing of patch data. Supports most easy types
    forAll(patchIDs, i)
    {
        label patchI = patchIDs[i];

        const globalIndex globalFaces(bf[patchI].size());

        if (isA<patchFieldType>(bf[patchI]))
        {
            // Explicit handling of externalCoupledObjectMixed bcs - they have
            // specialised writing routines

            const patchFieldType& pf = refCast<const patchFieldType>
            (
                bf[patchI]
            );
            OStringStream os;

            // Pass responsibility for all writing over to bc
            pf.writeData(os);

            // Collect contributions from all processors and output them on
            // master
            if (Pstream::master())
            {
                // Output master data first
                if (!headerDone)
                {
                    pf.writeHeader(masterFilePtr());
                    headerDone = true;
                }
                masterFilePtr() << os.str().c_str();

                for (label procI = 1; procI < Pstream::nProcs(); procI++)
                {
                    IPstream fromSlave(Pstream::scheduled, procI);
                    string str(fromSlave);
                    masterFilePtr() << str.c_str();
                }
            }
            else
            {
                OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
                toMaster << os.str();
            }
        }
        else if (isA<mixedFvPatchField<Type> >(bf[patchI]))
        {
            const mixedFvPatchField<Type>& pf =
                refCast<const mixedFvPatchField<Type> >(bf[patchI]);

            Field<Type> value(gatherAndCombine(pf));
            Field<Type> snGrad(gatherAndCombine(pf.snGrad()()));
            Field<Type> refValue(gatherAndCombine(pf.refValue()));
            Field<Type> refGrad(gatherAndCombine(pf.refGrad()));
            scalarField valueFraction(gatherAndCombine(pf.valueFraction()));

            if (Pstream::master())
            {
                forAll(refValue, faceI)
                {
                    masterFilePtr()
                        << value[faceI] << token::SPACE
                        << snGrad[faceI] << token::SPACE
                        << refValue[faceI] << token::SPACE
                        << refGrad[faceI] << token::SPACE
                        << valueFraction[faceI] << nl;
                }
            }
        }
        else
        {
            // Output the value and snGrad
            Field<Type> value(gatherAndCombine(bf[patchI]));
            Field<Type> snGrad(gatherAndCombine(bf[patchI].snGrad()()));
            if (Pstream::master())
            {
                forAll(value, faceI)
                {
                    masterFilePtr()
                        << value[faceI] << token::SPACE
                        << snGrad[faceI] << nl;
                }
            }
        }
    }

    return true;
}


// ************************************************************************* //
