/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "PatchFunction1.H"
#include "fieldTypes.H"
#include "ConstantField.H"
#include "UniformValueField.H"
#include "MappedFile.H"
#include "addToRunTimeSelectionTable.H"
#include "Table.H"

#define makePatchFunction1s(Type)                                              \
    makePatchFunction1(Type);                                                  \
    makePatchFunction1Type(ConstantField, Type);                               \
    makePatchFunction1Type(MappedFile, Type);                                  \
    makePatchFunction1Type(UniformValueField, Type);

#define addUniformValueFieldFunction1s(F1Type, Type)                           \
    PatchFunction1<Type>::adddictionaryConstructorToTable                      \
    <PatchFunction1Types::UniformValueField<Type>>                             \
        add##F1Type##UniformValueField##Type##ConstructorToTable_(#F1Type);

namespace Foam
{
    makePatchFunction1(label);
    makePatchFunction1Type(ConstantField, label);

    makePatchFunction1s(scalar);
    makePatchFunction1s(vector);
    makePatchFunction1s(sphericalTensor);
    makePatchFunction1s(symmTensor);
    makePatchFunction1s(tensor);


    //- Option1 : add UniformFieldValue under the same name as Function1
    //            See makeFunction1s.C. Note that we do not need
    //            Constant & Uniform
    addUniformValueFieldFunction1s(zero, scalar);
    addUniformValueFieldFunction1s(zero, vector);
    addUniformValueFieldFunction1s(zero, sphericalTensor);
    addUniformValueFieldFunction1s(zero, symmTensor);
    addUniformValueFieldFunction1s(zero, tensor);

    addUniformValueFieldFunction1s(one, scalar);
    addUniformValueFieldFunction1s(one, vector);
    addUniformValueFieldFunction1s(one, sphericalTensor);
    addUniformValueFieldFunction1s(one, symmTensor);
    addUniformValueFieldFunction1s(one, tensor);

    addUniformValueFieldFunction1s(polynomial, scalar);
    addUniformValueFieldFunction1s(polynomial, vector);
    addUniformValueFieldFunction1s(polynomial, sphericalTensor);
    addUniformValueFieldFunction1s(polynomial, symmTensor);
    addUniformValueFieldFunction1s(polynomial, tensor);

    addUniformValueFieldFunction1s(cosine, scalar);
    addUniformValueFieldFunction1s(cosine, vector);
    addUniformValueFieldFunction1s(cosine, sphericalTensor);
    addUniformValueFieldFunction1s(cosine, symmTensor);
    addUniformValueFieldFunction1s(cosine, tensor);

    addUniformValueFieldFunction1s(sine, scalar);
    addUniformValueFieldFunction1s(sine, vector);
    addUniformValueFieldFunction1s(sine, sphericalTensor);
    addUniformValueFieldFunction1s(sine, symmTensor);
    addUniformValueFieldFunction1s(sine, tensor);

    addUniformValueFieldFunction1s(square, scalar);
    addUniformValueFieldFunction1s(square, vector);
    addUniformValueFieldFunction1s(square, sphericalTensor);
    addUniformValueFieldFunction1s(square, symmTensor);
    addUniformValueFieldFunction1s(square, tensor);

    addUniformValueFieldFunction1s(csvFile, scalar);
    addUniformValueFieldFunction1s(csvFile, vector);
    addUniformValueFieldFunction1s(csvFile, sphericalTensor);
    addUniformValueFieldFunction1s(csvFile, symmTensor);
    addUniformValueFieldFunction1s(csvFile, tensor);

    addUniformValueFieldFunction1s(table, scalar);
    addUniformValueFieldFunction1s(table, vector);
    addUniformValueFieldFunction1s(table, sphericalTensor);
    addUniformValueFieldFunction1s(table, symmTensor);
    addUniformValueFieldFunction1s(table, tensor);

    addUniformValueFieldFunction1s(tableFile, scalar);
    addUniformValueFieldFunction1s(tableFile, vector);
    addUniformValueFieldFunction1s(tableFile, sphericalTensor);
    addUniformValueFieldFunction1s(tableFile, symmTensor);
    addUniformValueFieldFunction1s(tableFile, tensor);

    addUniformValueFieldFunction1s(scale, scalar);
    addUniformValueFieldFunction1s(scale, vector);
    addUniformValueFieldFunction1s(scale, sphericalTensor);
    addUniformValueFieldFunction1s(scale, symmTensor);
    addUniformValueFieldFunction1s(scale, tensor);


    ////- Option2 : at static initialisation add all Function1 types.
    ////  This does not work because we cannot guarantee that the Function1
    ////  static initialisation has happened already.
    //template<class Type>
    //class addToUniform
    //{
    //public:
    //    addToUniform()
    //    {
    //        // Get the Function1 table
    //        typedef typename Function1<Type>::dictionaryConstructorTable
    //            F1Type;
    //        Function1<Type>::constructdictionaryConstructorTables();
    //        const F1Type& F1Table =
    //            *Function1<Type>::dictionaryConstructorTablePtr_;
    //
    //        // Get the PatchFunction1 table
    //        typedef typename PatchFunction1<Type>::dictionaryConstructorTable
    //            PF1Type;
    //
    //        PatchFunction1<Type>::constructdictionaryConstructorTables();
    //        PF1Type& PF1Table =
    //            *PatchFunction1<Type>::dictionaryConstructorTablePtr_;
    //
    //        // Get the UniformValueField constructor
    //        auto cstrIter =
    //            PatchFunction1<Type>::dictionaryConstructorTablePtr_->cfind
    //            (
    //                PatchFunction1Types::UniformValueField<Type>::typeName
    //            );
    //
    //        // Add the UniformValueField under the Function1 name
    //        forAllConstIters(F1Table, iter)
    //        {
    //            //bool ok =
    //            PF1Table.insert(iter.key(), cstrIter());
    //            //if (!ok)
    //            //{
    //            //    std::cout<< "** problem" << std::endl;
    //            //}
    //        }
    //    }
    //};
    //static const addToUniform<scalar> addScalar;
    //static const addToUniform<vector> addVector;
    //static const addToUniform<sphericalTensor> addSphericalTensor;
    //static const addToUniform<symmTensor> addSymmTensor;
    //static const addToUniform<tensor> addTensor;
}


// ************************************************************************* //
