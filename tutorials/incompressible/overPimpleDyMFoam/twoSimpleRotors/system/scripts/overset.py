from paraview.simple import *
from paraview import coprocessing

#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# paraview version 5.5.0

#--------------------------------------------------------------
# Global screenshot output options
imageFileNamePadding=4
rescale_lookuptable=False

# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      # trace generated using paraview version 5.5.0

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [1077, 763]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [0.00784385809674859, 0.005000000004656613, 0.004999999888241291]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [0.0072242101003740155, 0.0002877833685303474, 0.035060283710920806]
      renderView1.CameraFocalPoint = [0.00868966107678934, 0.004150999005211765, 0.0049322758242629034]
      renderView1.CameraViewUp = [0.3542102656908786, 0.9252429122682538, 0.135869941401907]
      renderView1.CameraParallelScale = 0.00787069031419879
      renderView1.CameraParallelProjection = 1
      renderView1.Background = [0.32, 0.34, 0.43]

      # init the 'GridAxes3DActor' selected for 'AxesGrid'
      renderView1.AxesGrid.XTitleFontFile = ''
      renderView1.AxesGrid.YTitleFontFile = ''
      renderView1.AxesGrid.ZTitleFontFile = ''
      renderView1.AxesGrid.XLabelFontFile = ''
      renderView1.AxesGrid.YLabelFontFile = ''
      renderView1.AxesGrid.ZLabelFontFile = ''

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='insitu/image_%t.png', freq=1, fittoscreen=0, magnification=1, width=1077, height=763, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # restore active view
      SetActiveView(renderView1)
      # ----------------------------------------------------------------

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # a producer from a simulation input
      input1 = coprocessor.CreateProducer(datadescription, 'mesh')

      # cellMask [0,1]
      threshold1 = Threshold(Input=input1)
      threshold1.Scalars = ['CELLS', 'cellMask']
      threshold1.ThresholdRange = [0.9, 1.1]

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from threshold1
      threshold1Display = Show(threshold1, renderView1)

      # get color transfer function/color map for 'cellTypes'
      cellTypesLUT = GetColorTransferFunction('cellTypes')
      cellTypesLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 1.000244140625, 0.865003, 0.865003, 0.865003, 2.00048828125, 0.705882, 0.0156863, 0.14902]
      cellTypesLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'cellTypes'
      cellTypesPWF = GetOpacityTransferFunction('cellTypes')
      cellTypesPWF.Points = [0.0, 0.0, 0.5, 0.0, 2.00048828125, 1.0, 0.5, 0.0]
      cellTypesPWF.ScalarRangeInitialized = 1

      # trace defaults for the display properties.
      threshold1Display.Representation = 'Surface With Edges'
      threshold1Display.ColorArrayName = ['CELLS', 'cellTypes']
      threshold1Display.LookupTable = cellTypesLUT
      threshold1Display.OSPRayScaleArray = 'U'
      threshold1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      threshold1Display.SelectOrientationVectors = 'None'
      threshold1Display.ScaleFactor = 0.0019999999552965165
      threshold1Display.SelectScaleArray = 'None'
      threshold1Display.GlyphType = 'Arrow'
      threshold1Display.GlyphTableIndexArray = 'None'
      threshold1Display.GaussianRadius = 9.999999776482583e-05
      threshold1Display.SetScaleArray = ['POINTS', 'U']
      threshold1Display.ScaleTransferFunction = 'PiecewiseFunction'
      threshold1Display.OpacityArray = ['POINTS', 'U']
      threshold1Display.OpacityTransferFunction = 'PiecewiseFunction'
      threshold1Display.DataAxesGrid = 'GridAxesRepresentation'
      threshold1Display.SelectionCellLabelFontFile = ''
      threshold1Display.SelectionPointLabelFontFile = ''
      threshold1Display.PolarAxes = 'PolarAxesRepresentation'
      threshold1Display.ScalarOpacityFunction = cellTypesPWF
      threshold1Display.ScalarOpacityUnitDistance = 0.0017065741933059136

      # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
      threshold1Display.ScaleTransferFunction.Points = [-0.2505497634410858, 0.0, 0.5, 0.0, 0.3270378112792969, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
      threshold1Display.OpacityTransferFunction.Points = [-0.2505497634410858, 0.0, 0.5, 0.0, 0.3270378112792969, 1.0, 0.5, 0.0]

      # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
      threshold1Display.DataAxesGrid.XTitleFontFile = ''
      threshold1Display.DataAxesGrid.YTitleFontFile = ''
      threshold1Display.DataAxesGrid.ZTitleFontFile = ''
      threshold1Display.DataAxesGrid.XLabelFontFile = ''
      threshold1Display.DataAxesGrid.YLabelFontFile = ''
      threshold1Display.DataAxesGrid.ZLabelFontFile = ''

      # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
      threshold1Display.PolarAxes.PolarAxisTitleFontFile = ''
      threshold1Display.PolarAxes.PolarAxisLabelFontFile = ''
      threshold1Display.PolarAxes.LastRadialAxisTextFontFile = ''
      threshold1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for cellTypesLUT in view renderView1
      cellTypesLUTColorBar = GetScalarBar(cellTypesLUT, renderView1)
      cellTypesLUTColorBar.Title = 'cellTypes'
      cellTypesLUTColorBar.ComponentTitle = ''
      cellTypesLUTColorBar.TitleFontFile = ''
      cellTypesLUTColorBar.LabelFontFile = ''

      # set color bar visibility
      cellTypesLUTColorBar.Visibility = 1

      # show color legend
      threshold1Display.SetScalarBarVisibility(renderView1, True)

    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # Frequencies at which the coprocessor updates.
  freqs = {'mesh': [1, 1, 1]}
  coprocessor.SetUpdateFrequencies(freqs)
  return coprocessor


#--------------------------------------------------------------
# Global variable that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView and the update frequency
coprocessor.EnableLiveVisualization(True, 1)

# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor
    if datadescription.GetForceOutput() == True:
        # We are just going to request all fields and meshes from the simulation
        # code/adaptor.
        for i in range(datadescription.GetNumberOfInputDescriptions()):
            datadescription.GetInputDescription(i).AllFieldsOn()
            datadescription.GetInputDescription(i).GenerateMeshOn()
        return

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=rescale_lookuptable,
        image_quality=0, padding_amount=imageFileNamePadding)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
