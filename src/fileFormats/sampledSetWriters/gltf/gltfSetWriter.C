/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "gltfSetWriter.H"
#include "coordSet.H"
#include "fileName.H"
#include "OFstream.H"
#include "floatVector.H"
#include "foamGltfScene.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::Enum<typename Foam::gltfSetWriter<Type>::fieldOption>
Foam::gltfSetWriter<Type>::fieldOptionNames_
({
    { fieldOption::UNIFORM, "uniform" },
    { fieldOption::FIELD, "field" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::word Foam::gltfSetWriter<Type>::getColourMap
(
    const dictionary& dict
) const
{
    word colourMap = colourTable::predefinedNames.names()[0];
    dict.readIfPresent("colourMap", colourMap);

    return colourMap;
}


template<class Type>
const Foam::colourTable& Foam::gltfSetWriter<Type>::getColourTable
(
    const dictionary& dict
) const
{
    return colourTable::ref(getColourMap(dict));
}


template<class Type>
Foam::scalar Foam::gltfSetWriter<Type>::getFieldMin
(
    const word& fieldName
) const
{
    const dictionary fieldDict = fieldInfoDict_.subOrEmptyDict(fieldName);

    return fieldDict.getOrDefault("min", -GREAT);
}


template<class Type>
Foam::scalar Foam::gltfSetWriter<Type>::getFieldMax
(
    const word& fieldName
) const
{
    const dictionary fieldDict = fieldInfoDict_.subOrEmptyDict(fieldName);

    return fieldDict.getOrDefault("max", GREAT);
}


template<class Type>
Foam::tmp<Foam::scalarField> Foam::gltfSetWriter<Type>::getAlphaField
(
    const dictionary& dict,
    const wordList& valueSetNames,
    const List<const Field<Type>*>& valueSets
) const
{
    if (dict.found("alpha"))
    {
        const auto option = fieldOptionNames_.get("alpha", dict);

        switch (option)
        {
            case fieldOption::UNIFORM:
            {
                const scalar value = dict.getScalar("alphaValue");
                return tmp<scalarField>::New(valueSets[0]->size(), value);
            }
            case fieldOption::FIELD:
            {
                const word alphaFieldName = dict.get<word>("alphaField");
                const bool normalise = dict.get<bool>("normalise");
                const label i = valueSetNames.find(alphaFieldName);
                if (i == -1)
                {
                    FatalErrorInFunction
                        << "Unable to find field " << alphaFieldName
                        << ". Valid field names are:" << valueSetNames
                        << exit(FatalError);
                }

                auto tresult =
                    tmp<scalarField>::New(valueSets[i]->component(0));

                if (normalise)
                {
                    tresult.ref() /= mag(tresult() + ROOTVSMALL);
                }

                return tresult;
            }
        }
    }

    return tmp<scalarField>::New(valueSets[0]->size(), Zero);
}


template<class Type>
Foam::tmp<Foam::scalarField> Foam::gltfSetWriter<Type>::getTrackAlphaField
(
    const dictionary& dict,
    const wordList& valueSetNames,
    const List<List<Field<Type>>>& valueSets,
    const label tracki
) const
{
    if (dict.found("alpha"))
    {
        const auto option = fieldOptionNames_.get("alpha", dict);

        switch (option)
        {
            case fieldOption::UNIFORM:
            {
                const scalar value = dict.getScalar("alphaValue");
                return tmp<scalarField>::New
                (
                    valueSets[0][tracki].size(), value
                );
            }
            case fieldOption::FIELD:
            {
                const word alphaFieldName = dict.get<word>("alphaField");
                const bool normalise = dict.get<bool>("normalise");
                const label fieldi = valueSetNames.find(alphaFieldName);
                if (fieldi == -1)
                {
                    FatalErrorInFunction
                        << "Unable to find field " << alphaFieldName
                        << ". Valid field names are:" << valueSetNames
                        << exit(FatalError);
                }

                // Note: selecting the first component!
                auto tresult =
                    tmp<scalarField>::New
                    (
                        valueSets[fieldi][tracki].component(0)
                    );

                if (normalise)
                {
                    tresult.ref() /= mag(tresult() + ROOTVSMALL);
                }

                return tresult;
            }
        }
    }

    return tmp<scalarField>::New(valueSets[0][tracki].size(), Zero);
}


template<class Type>
Foam::vector Foam::gltfSetWriter<Type>::getTrackAnimationColour
(
    const colourTable& colours,
    const wordList& valueSetNames,
    const List<List<Field<Type>>>& valueSets,
    const label tracki
) const
{
    if (!colour_)
    {
        FatalErrorInFunction
            << "Attempting to get colour when colour option is off"
            << abort(FatalError);
    }

    const auto option = fieldOptionNames_.get("colour", animationDict_);

    switch (option)
    {
        case fieldOption::UNIFORM:
        {
            return animationDict_.get<vector>("colourValue");
        }
        case fieldOption::FIELD:
        {
            const word fieldName = animationDict_.get<word>("colourField");
            const label fieldi = valueSetNames.find(fieldName);
            if (fieldi == -1)
            {
                FatalErrorInFunction
                    << "Unable to find field " << fieldName
                    << ". Valid field names are:" << valueSetNames
                    << exit(FatalError);
            }

            // Note: selecting the first component!

            scalar minValue;
            scalar maxValue;
            if (!animationDict_.readIfPresent("min", minValue))
            {
                minValue = min(valueSets[fieldi][tracki].component(0));
            }
            if (!animationDict_.readIfPresent("max", maxValue))
            {
                maxValue = max(valueSets[fieldi][tracki].component(0));
            }
            const scalar refValue = component(valueSets[fieldi][tracki][0], 0);
            const scalar fraction =
                (refValue - minValue)/(maxValue - minValue + ROOTVSMALL);

            return (colours.value(max(0, min(1, fraction))));
        }
    }

    return vector::zero;
}


template<class Type>
Foam::tmp<Foam::vectorField> Foam::gltfSetWriter<Type>::directions
(
    const coordSet& points
) const
{
    auto tresult = tmp<vectorField>::New(points.size(), Zero);
    auto& result = tresult.ref();

    if (points.size() > 1)
    {
        for (label i = 1; i < points.size(); ++i)
        {
            result[i-1] = points[i] - points[i-1];
            result[i-1].normalise();
        }

        result.last() = result[points.size()-2];
    }


    return tresult;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::gltfSetWriter<Type>::gltfSetWriter()
:
    writer<Type>(),
    animate_(false),
    colour_(false),
    fieldInfoDict_(),
    animationDict_()
{}


template<class Type>
Foam::gltfSetWriter<Type>::gltfSetWriter(const dictionary& dict)
:
    writer<Type>(dict),
    animate_(dict.getOrDefault("animate", false)),
    colour_(dict.getOrDefault("colour", false)),
    fieldInfoDict_(dict.subOrEmptyDict("fieldInfo")),
    animationDict_(dict.subOrEmptyDict("animationInfo"))
{
    // fieldInfo
    // {
    //     U
    //     {
    //         colourMap       coolToWarm; // others...
    //         min             10;
    //         max             100;
    //         alpha           field; // uniform|field
    //         alphaField      ageOfAir;
    //     }
    // }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::gltfSetWriter<Type>::getFileName
(
    const coordSet& points,
    const wordList& valueSetNames
) const
{
    return this->getBaseName(points, valueSetNames) + ".gltf";
}


template<class Type>
void Foam::gltfSetWriter<Type>::write
(
    const coordSet& points,
    const wordList& valueSetNames,
    const List<const Field<Type>*>& valueSets,
    Ostream& os
) const
{
    if (valueSets.size() != valueSetNames.size())
    {
        FatalErrorInFunction
            << "Number of variables:" << valueSetNames.size() << endl
            << "Number of valueSets:" << valueSets.size()
            << exit(FatalError);
    }

    glTF::scene scene;
    const label meshi = scene.addMesh(points, "points");
    forAll(valueSetNames, i)
    {
        scene.addFieldToMesh(*valueSets[i], valueSetNames[i], meshi);
    }

    if (colour_)
    {
        forAll(valueSets, fieldi)
        {
            const auto& field = *valueSets[fieldi];
            const word& fieldName = valueSetNames[fieldi];
            const dictionary dict = fieldInfoDict_.subOrEmptyDict(fieldName);
            const auto& colours = getColourTable(dict);

            const auto talpha =
                 getAlphaField(dict, valueSetNames, valueSets);
            const scalarField& alpha = talpha();

            const Type maxValue = max(field);
            const Type minValue = min(field);

            const scalar minValueLimit = getFieldMin(fieldName);
            const scalar maxValueLimit = getFieldMax(fieldName);

            for (direction cmpti=0; cmpti < pTraits<Type>::nComponents; ++cmpti)
            {
                vectorField fieldColour(field.size());

                forAll(field, i)
                {
                    const Type& v = field[i];
                    float f = component(v, cmpti);
                    float minf = max(component(minValue, cmpti), minValueLimit);
                    float maxf = min(component(maxValue, cmpti), maxValueLimit);
                    float deltaf = (maxf - minf + SMALL);

                    fieldColour[i] =
                        colours.value(min(max((f - minf)/deltaf, 0), 1));
                }

                scene.addColourToMesh
                (
                    fieldColour,
                    "Colour:" + fieldName + Foam::name(cmpti),
                    meshi,
                    alpha
                );
            }
        }
    }

    scene.write(os);
}


template<class Type>
void Foam::gltfSetWriter<Type>::write
(
    const bool writeTracks,
    const List<scalarField>& times,
    const PtrList<coordSet>& tracks,
    const wordList& valueSetNames,
    const List<List<Field<Type>>>& valueSets,
    Ostream& os
) const
{
    if (valueSets.size() != valueSetNames.size())
    {
        FatalErrorInFunction
            << "Number of variables:" << valueSetNames.size() << endl
            << "Number of valueSets:" << valueSets.size()
            << exit(FatalError);
    }

    if (animate_)
    {
        writeAnimateTracks
        (
            writeTracks,
            times,
            tracks,
            valueSetNames,
            valueSets,
            os
        );
    }
    else
    {
        writeStaticTracks
        (
            writeTracks,
            times,
            tracks,
            valueSetNames,
            valueSets,
            os
        );
    }
}


template<class Type>
void Foam::gltfSetWriter<Type>::writeStaticTracks
(
    const bool writeTracks,
    const List<scalarField>& times,
    const PtrList<coordSet>& tracks,
    const wordList& valueSetNames,
    const List<List<Field<Type>>>& valueSets,
    Ostream& os
) const
{
    glTF::scene scene;
    forAll(tracks, tracki)
    {
        const vectorField& track = tracks[tracki];
        const label meshi = scene.addMesh(track, "track:" + Foam::name(tracki));
        forAll(valueSetNames, fieldi)
        {
            const word& fieldName = valueSetNames[fieldi];
            const auto& field = valueSets[fieldi][tracki];
            scene.addFieldToMesh(field, fieldName, meshi);
        }

        if (colour_)
        {
            forAll(valueSets, fieldi)
            {
                const auto& field = valueSets[fieldi][tracki];
                const word& fieldName = valueSetNames[fieldi];
                const dictionary dict =
                    fieldInfoDict_.subOrEmptyDict(fieldName);
                const auto& colours = getColourTable(dict);

                const auto talpha =
                    getTrackAlphaField(dict, valueSetNames, valueSets, tracki);
                const scalarField& alpha = talpha();

                const Type maxValue = max(field);
                const Type minValue = min(field);

                const scalar minValueLimit = getFieldMin(fieldName);
                const scalar maxValueLimit = getFieldMax(fieldName);

                for
                (
                    direction cmpti=0;
                    cmpti < pTraits<Type>::nComponents;
                    ++cmpti
                )
                {
                    vectorField fieldColour(field.size(), Zero);

                    forAll(field, i)
                    {
                        const Type& v = field[i];
                        float f = component(v, cmpti);
                        float minf =
                            max(component(minValue, cmpti), minValueLimit);
                        float maxf =
                            min(component(maxValue, cmpti), maxValueLimit);
                        float deltaf = (maxf - minf + SMALL);

                        fieldColour[i] =
                           colours.value(min(max((f - minf)/deltaf, 0), 1));
                    }

                    scene.addColourToMesh
                    (
                        fieldColour,
                        "Colour:" + fieldName + Foam::name(cmpti),
                        meshi,
                        alpha
                    );
                }
            }
        }
    }

    scene.write(os);
}


template<class Type>
void Foam::gltfSetWriter<Type>::writeAnimateTracks
(
    const bool writeTracks,
    const List<scalarField>& times,
    const PtrList<coordSet>& tracks,
    const wordList& valueSetNames,
    const List<List<Field<Type>>>& valueSets,
    Ostream& os
) const
{
    const auto& colours = getColourTable(animationDict_);

    glTF::scene scene;
    const label animationi = scene.createAnimation("animation");

    forAll(tracks, tracki)
    {
        const auto& track = tracks[tracki];

        if (track.empty())
        {
            continue;
        }

        // Seed starting positions and field values
        const label meshi =
            scene.addMesh
            (
                vectorField(1, track[0]),
                "track:" + Foam::name(tracki)
            );

        forAll(valueSetNames, fieldi)
        {
            const Field<Type>& field = valueSets[fieldi][tracki];
            const word& fieldName = valueSetNames[fieldi];
            scene.addFieldToMesh(Field<Type>(1, field[0]), fieldName, meshi);
        }

        // Time frames
        const label timeId =
            scene.addField(times[tracki], "time:" + Foam::name(tracki));

        // Translations
        const vectorField translation(track - track[0]);
        const label translationId = scene.addField(translation, "translation");

        scene.addToAnimation(animationi, timeId, translationId, meshi);

        // Note: colours cannot be animated... setting a fixed value
        if (colour_)
        {
            const vector colour =
                getTrackAnimationColour
                (
                    colours,
                    valueSetNames,
                    valueSets,
                    tracki
                );

            const auto talpha =
                getTrackAlphaField
                (
                    animationDict_,
                    valueSetNames,
                    valueSets,
                    tracki
                );

            const scalarField& alpha = talpha();

            scene.addColourToMesh
            (
                vectorField(1, colour),
                "Colour:fixed",
                meshi,
                scalarField(1, alpha[0])
            );
        }
    }

    scene.write(os);
}


// ************************************************************************* //
