/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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

// OpenFOAM includes
#include "fieldVisualisationBase.H"
#include "runTimePostProcessing.H"

// VTK includes
#include "vtkArrowSource.h"
#include "vtkCellData.h"
#include "vtkColorTransferFunction.h"
#include "vtkFloatArray.h"
#include "vtkGlyph3D.h"
#include "vtkLookupTable.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderer.h"
#include "vtkScalarBarActor.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkCellDataToPointData.h"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::functionObjects::fieldVisualisationBase::colourByType
>
Foam::functionObjects::fieldVisualisationBase::colourByTypeNames
{
    { colourByType::cbColour, "colour" },
    { colourByType::cbField, "field" },
};

const Foam::Enum
<
    Foam::functionObjects::fieldVisualisationBase::colourMapType
>
Foam::functionObjects::fieldVisualisationBase::colourMapTypeNames
{
    { colourMapType::cmRainbow, "rainbow" },
    { colourMapType::cmBlueWhiteRed, "blueWhiteRed" },
    { colourMapType::cmFire, "fire" },
    { colourMapType::cmGreyscale, "greyscale" },
};


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fieldVisualisationBase::setColourMap
(
    vtkLookupTable* lut
) const
{
    label nColours = 256;

    lut->SetNumberOfColors(nColours);

    auto ctf = vtkSmartPointer<vtkColorTransferFunction>::New();

    switch (colourMap_)
    {
        case cmRainbow:
        {
            ctf->SetColorSpaceToHSV();
            ctf->AddRGBPoint(0, 0, 0, 1);
            ctf->AddRGBPoint(0.5, 0, 1, 0);
            ctf->AddRGBPoint(1, 1, 0, 0);
            break;
        }
        case cmBlueWhiteRed:
        {
            // Values taken from ParaView settings
            ctf->SetColorSpaceToDiverging();
            ctf->AddRGBPoint(0.0, 0.231373, 0.298039, 0.752941);
            ctf->AddRGBPoint(0.5, 0.865003, 0.865003, 0.865003);
            ctf->AddRGBPoint(1.0, 0.705882, 0.0156863, 0.14902);
            break;
        }
        case cmFire:
        {
            // Values taken from ParaView settings
            ctf->SetColorSpaceToRGB();
            ctf->AddRGBPoint(0, 0, 0, 0);
            ctf->AddRGBPoint(0.4, 0.901961, 0, 0);
            ctf->AddRGBPoint(0.8, 0.901961, 0.901961, 0);
            ctf->AddRGBPoint(1, 1, 1, 1);
            break;
        }
        case cmGreyscale:
        {
            ctf->SetColorSpaceToRGB();
            ctf->AddRGBPoint(0, 0, 0, 0);
            ctf->AddRGBPoint(1, 1, 1, 1);
            break;
        }
    }


    for (label i = 0; i < nColours; i++)
    {
        double* c = ctf->GetColor(scalar(i)/scalar(nColours));
        lut->SetTableValue(i, c[0], c[1], c[2], 1.0);
    }
}


void Foam::functionObjects::fieldVisualisationBase::addScalarBar
(
    const scalar position,
    vtkRenderer* renderer,
    vtkLookupTable* lut
) const
{
    // Add scalar bar legend
    if (!scalarBar_.visible_)
    {
        return;
    }

    auto sbar = vtkSmartPointer<vtkScalarBarActor>::New();
    sbar->SetLookupTable(lut);
    sbar->SetNumberOfLabels(scalarBar_.numberOfLabels_);

    const vector textColour = colours_["text"]->value(position);

    // Work-around to supply our own scalarbar title
    // - Default scalar bar title text is scales by the scalar bar box
    //   dimensions so if the title is a long string, the text is shrunk to fit
    //   Instead, suppress title and set the title using a vtkTextActor
    auto titleActor = vtkSmartPointer<vtkTextActor>::New();
    sbar->SetTitle(" ");
    titleActor->SetInput(scalarBar_.title_.c_str());
    titleActor->GetTextProperty()->SetFontFamilyToArial();
    titleActor->GetTextProperty()->SetFontSize(3*scalarBar_.fontSize_);
    titleActor->GetTextProperty()->SetJustificationToCentered();
    titleActor->GetTextProperty()->SetVerticalJustificationToBottom();
    titleActor->GetTextProperty()->BoldOn();
    titleActor->GetTextProperty()->ItalicOff();
    titleActor->GetTextProperty()->SetColor
    (
        textColour[0],
        textColour[1],
        textColour[2]
    );
    titleActor->GetPositionCoordinate()->
        SetCoordinateSystemToNormalizedViewport();

    // How to use the standard scalar bar text
    // sbar->SetTitle(scalarBar_.title_.c_str());
    // sbar->GetTitleTextProperty()->SetColor
    // (
    //     textColour[0],
    //     textColour[1],
    //     textColour[2]
    // );
    // sbar->GetTitleTextProperty()->SetFontSize(scalarBar_.fontSize_);
    // sbar->GetTitleTextProperty()->ShadowOff();
    // sbar->GetTitleTextProperty()->BoldOn();
    // sbar->GetTitleTextProperty()->ItalicOff();

    sbar->GetLabelTextProperty()->SetColor
    (
        textColour[0],
        textColour[1],
        textColour[2]
    );
    sbar->GetLabelTextProperty()->SetFontSize(scalarBar_.fontSize_);
    sbar->GetLabelTextProperty()->ShadowOff();
    sbar->GetLabelTextProperty()->BoldOff();
    sbar->GetLabelTextProperty()->ItalicOff();
    sbar->SetLabelFormat(scalarBar_.labelFormat_.c_str());

    sbar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    sbar->GetPositionCoordinate()->SetValue
    (
        scalarBar_.position_.first(),
        scalarBar_.position_.second()
    );
    if (scalarBar_.vertical_)
    {
        sbar->SetOrientationToVertical();
        sbar->SetWidth(0.1);
        sbar->SetHeight(0.75);
        sbar->SetTextPositionToSucceedScalarBar();
    }
    else
    {
        sbar->SetOrientationToHorizontal();

        // Adjustments since not using scalarbar title property
        sbar->SetWidth(0.75);
        sbar->SetHeight(0.07);
        sbar->SetBarRatio(0.5);
        // sbar->SetHeight(0.1);
        // sbar->SetTitleRatio(0.01);
        sbar->SetTextPositionToPrecedeScalarBar();
    }

    titleActor->GetPositionCoordinate()->SetValue
    (
        scalarBar_.position_.first() + 0.5*sbar->GetWidth(),
        scalarBar_.position_.second() + sbar->GetHeight()
    );

    // sbar->DrawFrameOn();
    // sbar->DrawBackgroundOn();
    // sbar->UseOpacityOff();
    // sbar->VisibilityOff();
    sbar->VisibilityOn();

    renderer->AddActor(sbar);
    renderer->AddActor2D(titleActor);
}


void Foam::functionObjects::fieldVisualisationBase::setField
(
    const scalar position,
    const word& colourFieldName,
    vtkPolyDataMapper* mapper,
    vtkRenderer* renderer,
    vtkPolyData* pData
) const
{
    mapper->InterpolateScalarsBeforeMappingOn();

    switch (colourBy_)
    {
        case cbColour:
        {
            mapper->ScalarVisibilityOff();
            break;
        }
        case cbField:
        {
            // Create look-up table for colours
            auto lut = vtkSmartPointer<vtkLookupTable>::New();
            setColourMap(lut);
            lut->SetVectorMode(vtkScalarsToColors::MAGNITUDE);

            // Configure the mapper
            const char* fieldName = colourFieldName.c_str();
            mapper->SelectColorArray(fieldName);

            // Set to use either point or cell data
            // Note: if both point and cell data exists, preferentially
            //       choosing point data.  This is often the case when using
            //       glyphs
            if (pData->GetPointData()->HasArray(fieldName) == 1)
            {
                mapper->SetScalarModeToUsePointFieldData();
            }
            else if (pData->GetCellData()->HasArray(fieldName) == 1)
            {
                mapper->SetScalarModeToUseCellFieldData();
            }
            else
            {
                WarningInFunction
                    << "Unable to determine cell or point data type "
                    << "- assuming point data";
                mapper->SetScalarModeToUsePointFieldData();
            }
            mapper->SetScalarRange(range_.first(), range_.second());
            mapper->SetColorModeToMapScalars();
            mapper->SetLookupTable(lut);
            mapper->ScalarVisibilityOn();

            // Add the bar
            addScalarBar(position, renderer, lut);
            break;
        }
    }

    mapper->Modified();
}



void Foam::functionObjects::fieldVisualisationBase::addGlyphs
(
    const scalar position,
    const word& scaleFieldName,
    const word& colourFieldName,
    const scalar maxGlyphLength,
    vtkPolyData* data,
    vtkActor* actor,
    vtkRenderer* renderer
) const
{
    auto glyph = vtkSmartPointer<vtkGlyph3D>::New();
    auto glyphMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    glyphMapper->SetInputConnection(glyph->GetOutputPort());

    glyph->SetInputData(data);
    glyph->ScalingOn();
    bool ok = true;

    // Determine whether we have scalar or vector data
    label nComponents = -1;
    const char* scaleFieldNameChar = scaleFieldName.c_str();
    if (data->GetPointData()->HasArray(scaleFieldNameChar) == 1)
    {
        nComponents =
            data->GetPointData()->GetArray(scaleFieldNameChar)
                ->GetNumberOfComponents();
    }
    else if (data->GetCellData()->HasArray(scaleFieldNameChar) == 1)
    {
        // Need to convert cell data to point data
        auto cellToPoint = vtkSmartPointer<vtkCellDataToPointData>::New();
        cellToPoint->SetInputData(data);
        cellToPoint->Update();
        vtkDataSet* pds = cellToPoint->GetOutput();
        vtkDataArray* pData = pds->GetPointData()->GetArray(scaleFieldNameChar);

        // Store in main vtkPolyData
        data->GetPointData()->AddArray(pData);

        nComponents = pData->GetNumberOfComponents();
    }
    else
    {
        WarningInFunction
            << "Glyphs can only be added to scalar or vector data. "
            << "Unable to process field " << scaleFieldName << endl;
        return;
    }

    if (nComponents == 1)
    {
        auto sphere = vtkSmartPointer<vtkSphereSource>::New();
        sphere->SetCenter(0, 0, 0);
        sphere->SetRadius(0.5);

        // Setting higher resolution slows the rendering significantly
        // sphere->SetPhiResolution(20);
        // sphere->SetThetaResolution(20);

        glyph->SetSourceConnection(sphere->GetOutputPort());

        if (maxGlyphLength > 0)
        {
            double range[2];

            // Can use values to find range
            // vtkDataArray* values =
            //     data->GetPointData()->GetScalars(scaleFieldNameChar);
            // values->GetRange(range);

            // Set range accoding to user-supplied limits
            range[0] = range_.first();
            range[1] = range_.second();
            glyph->ClampingOn();
            glyph->SetRange(range);

            // If range[0] != min(value), maxGlyphLength behaviour will not
            // be correct...
            glyph->SetScaleFactor(maxGlyphLength);
        }
        else
        {
            glyph->SetScaleFactor(1);
        }
        glyph->SetScaleModeToScaleByScalar();
        glyph->OrientOff();
        glyph->SetColorModeToColorByScalar();
        glyph->SetInputArrayToProcess
        (
            0, // scalars
            0,
            0,
            vtkDataObject::FIELD_ASSOCIATION_POINTS,
            scaleFieldNameChar
        );
    }
    else if (nComponents == 3)
    {
        auto arrow = vtkSmartPointer<vtkArrowSource>::New();
        arrow->SetTipResolution(10);
        arrow->SetTipRadius(0.1);
        arrow->SetTipLength(0.35);
        arrow->SetShaftResolution(10);
        arrow->SetShaftRadius(0.03);

        glyph->SetSourceConnection(arrow->GetOutputPort());

        if (maxGlyphLength > 0)
        {
            vtkDataArray* values =
                data->GetPointData()->GetVectors(scaleFieldNameChar);

            double range[6];
            values->GetRange(range);

            // Attempt to set range for vectors...
            // scalar x0 = sqrt(sqr(range_.first())/3.0);
            // scalar x1 = sqrt(sqr(range_.second())/3.0);
            // range[0] = x0;
            // range[1] = x0;
            // range[2] = x0;
            // range[3] = x1;
            // range[4] = x1;
            // range[5] = x1;

            glyph->ClampingOn();
            glyph->SetRange(range);
            glyph->SetScaleFactor(maxGlyphLength);
        }
        else
        {
            glyph->SetScaleFactor(1);
        }
        glyph->SetScaleModeToScaleByVector();
        glyph->OrientOn();
        glyph->SetVectorModeToUseVector();
        glyph->SetColorModeToColorByVector();
        glyph->SetInputArrayToProcess
        (
            1, // vectors
            0,
            0,
            vtkDataObject::FIELD_ASSOCIATION_POINTS,
            scaleFieldNameChar
        );
    }
    else
    {
        WarningInFunction
            << "Glyphs can only be added to scalar and vector fields."
            << " Field " << scaleFieldName << " has " << nComponents
            << " components" << endl;

        ok = false;
    }

    if (ok)
    {
        glyph->Update();

        setField(position, colourFieldName, glyphMapper, renderer, data);

        glyphMapper->Update();

        actor->SetMapper(glyphMapper);

        renderer->AddActor(actor);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldVisualisationBase::fieldVisualisationBase
(
    const dictionary& dict,
    const HashPtrTable<Function1<vector>, word>& colours
)
:
    colours_(colours),
    fieldName_(dict.lookup("field")),
    colourBy_(cbColour),
    colourMap_(cmRainbow),
    range_()
{
    colourBy_ = colourByTypeNames.lookup("colourBy", dict);

    switch (colourBy_)
    {
        case cbColour:
        {
            scalarBar_.visible_ = false;
            break;
        }
        case cbField:
        {
            dict.lookup("range") >> range_;

            if (dict.found("colourMap"))
            {
                colourMap_ = colourMapTypeNames.lookup("colourMap", dict);
            }

            const dictionary& sbarDict = dict.subDict("scalarBar");
            sbarDict.lookup("visible") >> scalarBar_.visible_;
            if (scalarBar_.visible_)
            {
                sbarDict.lookup("vertical") >> scalarBar_.vertical_;
                sbarDict.lookup("position") >> scalarBar_.position_;
                sbarDict.lookup("title") >> scalarBar_.title_;
                sbarDict.lookup("fontSize") >> scalarBar_.fontSize_;
                sbarDict.lookup("labelFormat") >> scalarBar_.labelFormat_;
                sbarDict.lookup("numberOfLabels") >> scalarBar_.numberOfLabels_;
            }
            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldVisualisationBase::~fieldVisualisationBase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::HashPtrTable<Foam::Function1<Foam::vector>, Foam::word>&
Foam::functionObjects::fieldVisualisationBase::colours() const
{
    return colours_;
}


const Foam::word&
Foam::functionObjects::fieldVisualisationBase::fieldName() const
{
    return fieldName_;
}


// ************************************************************************* //
