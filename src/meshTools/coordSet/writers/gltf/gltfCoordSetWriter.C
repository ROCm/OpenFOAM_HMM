/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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
    // No naming for NONE
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
Foam::scalarMinMax Foam::gltfSetWriter<Type>::getFieldLimits
(
    const word& fieldName
) const
{
    const dictionary fieldDict = fieldInfoDict_.subOrEmptyDict(fieldName);

    scalarMinMax limits;

    fieldDict.readIfPresent("min", limits.min());
    fieldDict.readIfPresent("max", limits.max());

    return limits;
}


template<class Type>
Foam::tmp<Foam::scalarField> Foam::gltfSetWriter<Type>::getAlphaField
(
    const dictionary& dict
) const
{
    // Fallback value
    scalar alphaValue(1);

    const entry* eptr = dict.findEntry("alpha", keyType::LITERAL);

    if (!eptr)
    {
        // Not specified
    }
    else if (!eptr->stream().peek().isString())
    {
        // Value specified

        ITstream& is = eptr->stream();
        is >> alphaValue;
        dict.checkITstream(is, "alpha");
    }
    else
    {
        // Enumeration

        const auto option = fieldOptionNames_.get("alpha", dict);

        switch (option)
        {
            case fieldOption::NONE:
            {
                break;
            }
            case fieldOption::UNIFORM:
            {
                dict.readEntry("alphaValue", alphaValue);
                break;
            }
            case fieldOption::FIELD:
            {
                WarningInFunction
                    << "Unsupported 'field' specification for alpha values"
                    << endl;
                break;
            }
        }
    }

    return tmp<scalarField>::New(1, alphaValue);
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

    const dictionary& dict = animationDict_;

    // Fallback value
    vector colourValue(Zero);

    const entry* eptr = dict.findEntry("colour", keyType::LITERAL);

    if (!eptr || !eptr->isStream())
    {
        FatalIOErrorInFunction(dict)
            << "Missing 'colour' entry"
            << exit(FatalIOError);
    }
    else if (!eptr->stream().peek().isString())
    {
        // Value specified

        ITstream& is = eptr->stream();
        is >> colourValue;
        dict.checkITstream(is, "colour");
    }
    else
    {
        // Enumeration

        const auto option = fieldOptionNames_.get("colour", dict);

        switch (option)
        {
            case fieldOption::NONE:
            {
                break;
            }
            case fieldOption::UNIFORM:
            {
                dict.readEntry("colourValue", colourValue);
                break;
            }
            case fieldOption::FIELD:
            {
                const word fieldName = dict.get<word>("colourField");
                const label fieldi = valueSetNames.find(fieldName);
                if (fieldi == -1)
                {
                    FatalErrorInFunction
                        << "Unable to find field " << fieldName
                        << ". Valid field names are:" << valueSetNames
                        << exit(FatalError);
                }

                const Field<Type>& colourFld = valueSets[fieldi][tracki];


                scalar refValue(0);
                scalarMinMax valLimits;

                if (pTraits<Type>::nComponents == 1)
                {
                    MinMax<Type> scanned(minMax(colourFld));

                    refValue = scalar(component(colourFld[0], 0));
                    valLimits.min() = scalar(component(scanned.min(), 0));
                    valLimits.max() = scalar(component(scanned.max(), 0));
                }
                else
                {
                    // Use mag() for multiple components
                    refValue = mag(colourFld[0]);
                    valLimits = minMaxMag(colourFld);
                }

                dict.readIfPresent("min", valLimits.min());
                dict.readIfPresent("max", valLimits.max());

                const scalar fraction =
                (
                    (refValue - valLimits.min())
                  / (valLimits.max() - valLimits.min() + ROOTVSMALL)
                );

                return colours.value(fraction);  // 0-1 clipped by value()
            }
        }
    }

    return colourValue;
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
    //         alpha           0.5;
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

            const auto talpha = getAlphaField(dict);
            const scalarField& alpha = talpha();

            const scalarMinMax valLimits = getFieldLimits(fieldName);

            // Generated field colours
            vectorField fieldColour(field.size());

            scalarMinMax fldLimits;

            if (pTraits<Type>::nComponents == 1)
            {
                MinMax<Type> scanned(minMax(field));

                fldLimits.min() = scalar(component(scanned.min(), 0));
                fldLimits.max() = scalar(component(scanned.max(), 0));
            }
            else
            {
                // Use mag() for multiple components
                fldLimits = minMaxMag(field);
            }

            const scalar minf = max(fldLimits.min(), valLimits.min());
            const scalar maxf = min(fldLimits.max(), valLimits.max());
            const scalar deltaf = (maxf - minf + SMALL);

            forAll(field, i)
            {
                const Type& val = field[i];

                const scalar f =
                (
                    pTraits<Type>::nComponents == 1
                  ? scalar(component(val, 0))
                  : scalar(mag(val))
                );

                // 0-1 clipped by value()
                fieldColour[i] = colours.value((f - minf)/deltaf);
            }

            scene.addColourToMesh
            (
                fieldColour,
                "Colour:" + fieldName,
                meshi,
                alpha
            );
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

                const auto talpha = getAlphaField(dict);
                const scalarField& alpha = talpha();

                const scalarMinMax valLimits = getFieldLimits(fieldName);


                // Generated field colours
                vectorField fieldColour(field.size());

                scalarMinMax fldLimits;

                if (pTraits<Type>::nComponents == 1)
                {
                    MinMax<Type> scanned(minMax(field));

                    fldLimits.min() = scalar(component(scanned.min(), 0));
                    fldLimits.max() = scalar(component(scanned.max(), 0));
                }
                else
                {
                    // Use mag() for multiple components
                    fldLimits = minMaxMag(field);
                }

                const scalar minf = max(fldLimits.min(), valLimits.min());
                const scalar maxf = min(fldLimits.max(), valLimits.max());
                const scalar deltaf = (maxf - minf + SMALL);

                forAll(field, i)
                {
                    const Type& val = field[i];

                    const scalar f =
                    (
                        pTraits<Type>::nComponents == 1
                      ? scalar(component(val, 0))
                      : scalar(mag(val))
                    );

                    // 0-1 clipped by value()
                    fieldColour[i] = colours.value((f - minf)/deltaf);
                }

                scene.addColourToMesh
                (
                    fieldColour,
                    "Colour:" + fieldName,
                    meshi,
                    alpha
                );
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

            const auto talpha = getAlphaField(animationDict_);

            const scalarField& alpha = talpha();

            scene.addColourToMesh
            (
                vectorField(1, colour),
                "Colour:fixed",  // ... or "Colour:constant"
                meshi,
                scalarField(1, alpha[0])
            );
        }
    }

    scene.write(os);
}


// ************************************************************************* //
