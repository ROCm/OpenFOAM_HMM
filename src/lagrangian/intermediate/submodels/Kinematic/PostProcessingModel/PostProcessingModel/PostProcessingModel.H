/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

Class
    Foam::PostProcessingModel

Description
    Templated post-processing model class

SourceFiles
    PostProcessingModel.C
    PostProcessingModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef PostProcessingModel_H
#define PostProcessingModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "SubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class PostProcessingModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PostProcessingModel
:
    public SubModelBase<CloudType>
{
    // Protected Member Functions

        //- Write post-processing info
        virtual void write();


public:

    //- Runtime type information
    TypeName("PostProcessingModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        PostProcessingModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        PostProcessingModel(CloudType& owner);

        //- Construct from dictionary
        PostProcessingModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        PostProcessingModel(const PostProcessingModel<CloudType>& ppm);

        //- Construct and return a clone
        virtual autoPtr<PostProcessingModel<CloudType> > clone() const
        {
            return autoPtr<PostProcessingModel<CloudType> >
            (
                new PostProcessingModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~PostProcessingModel();


    //- Selector
    static autoPtr<PostProcessingModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        // Evaluation

            //- Main post-processing function
            virtual void post();

            //- Gather post-processing info on patch
            virtual void postPatch
            (
                const typename CloudType::parcelType& p,
                const label patchI
            );

            //- Gather post-processing info on a face
            virtual void postFace(const typename CloudType::parcelType& p);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makePostProcessingModel(CloudType)                                    \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(PostProcessingModel<CloudType>, 0);   \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        PostProcessingModel<CloudType>,                                       \
        dictionary                                                            \
    );


#define makePostProcessingModelType(SS, CloudType, ParcelType)                \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType<ParcelType> >, 0);       \
                                                                              \
    PostProcessingModel<CloudType<ParcelType> >::                             \
        adddictionaryConstructorToTable<SS<CloudType<ParcelType> > >          \
            add##SS##CloudType##ParcelType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PostProcessingModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
