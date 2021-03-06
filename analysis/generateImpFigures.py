#-----------------------------------------------------------
# Helper script to generate FTLE figures for the rigid clots
#
# Author:       Zachariah Irwin
# Institution:  University of Colroado Boulder
# Last Edits:   May 2020
#-----------------------------------------------------------
import sys, os
from paraview.simple import *

structNames = ['M1', 'M2', 'M3', 'M4', 'M7', 'M8', 'M9']
leakNames   = ['0PL', '10PL', '20PL', '40PL']
timeIndexes = ['6','11','14','24','30','38','50','68']
timeNames   = ['T1','T2','T3','T4','T5','T6','T7','T8']

baseFTLEDir = os.getcwd() + '/output/'

for timeInd in range(len(timeIndexes)):

      ftleFileName = baseFTLEDir + 'FTLE-HS-Impermeable/FTLE-HS-Impermeable-%s.vtk' %(timeIndexes[timeInd])

      outName      = baseFTLEDir + 'FTLE_IMP/High-Shear-Clot/FTLE-HS-Impermeable-%s-.png' %(timeNames[timeInd])
      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      LoadPalette(paletteName='WhiteBackground')

      # get active view
      renderView1 = GetActiveViewOrCreate('RenderView')
      # uncomment following to set a specific view size
      # renderView1.ViewSize = [1483, 708]
      renderView1.UseLight = 0

      # Properties modified on renderView1
      renderView1.OrientationAxesVisibility = 0

      # get the material library
      materialLibrary1 = GetMaterialLibrary()

      # create a new 'Legacy VTK Reader'
      FTLE_FILE_vtk = LegacyVTKReader(FileNames=[ftleFileName])

      # show data in view
      FTLE_FILE_vtkDisplay = Show(FTLE_FILE_vtk, renderView1)

      # trace defaults for the display properties.
      FTLE_FILE_vtkDisplay.Representation = 'Surface'
      FTLE_FILE_vtkDisplay.AmbientColor = [0.0, 0.0, 0.0]
      FTLE_FILE_vtkDisplay.ColorArrayName = [None, '']
      FTLE_FILE_vtkDisplay.DiffuseColor = [1.0, 1.0, 1.0]
      FTLE_FILE_vtkDisplay.LookupTable = None
      FTLE_FILE_vtkDisplay.MapScalars = 1
      FTLE_FILE_vtkDisplay.MultiComponentsMapping = 0
      FTLE_FILE_vtkDisplay.InterpolateScalarsBeforeMapping = 1
      FTLE_FILE_vtkDisplay.Opacity = 1.0
      FTLE_FILE_vtkDisplay.PointSize = 2.0
      FTLE_FILE_vtkDisplay.LineWidth = 1.0
      FTLE_FILE_vtkDisplay.RenderLinesAsTubes = 0
      FTLE_FILE_vtkDisplay.RenderPointsAsSpheres = 0
      FTLE_FILE_vtkDisplay.Interpolation = 'Gouraud'
      FTLE_FILE_vtkDisplay.Specular = 0.0
      FTLE_FILE_vtkDisplay.SpecularColor = [1.0, 1.0, 1.0]
      FTLE_FILE_vtkDisplay.SpecularPower = 100.0
      FTLE_FILE_vtkDisplay.Luminosity = 0.0
      FTLE_FILE_vtkDisplay.Ambient = 0.0
      FTLE_FILE_vtkDisplay.Diffuse = 1.0
      FTLE_FILE_vtkDisplay.EdgeColor = [0.0, 0.0, 0.5]
      FTLE_FILE_vtkDisplay.BackfaceRepresentation = 'Follow Frontface'
      FTLE_FILE_vtkDisplay.BackfaceAmbientColor = [1.0, 1.0, 1.0]
      FTLE_FILE_vtkDisplay.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
      FTLE_FILE_vtkDisplay.BackfaceOpacity = 1.0
      FTLE_FILE_vtkDisplay.Position = [0.0, 0.0, 0.0]
      FTLE_FILE_vtkDisplay.Scale = [1.0, 1.0, 1.0]
      FTLE_FILE_vtkDisplay.Orientation = [0.0, 0.0, 0.0]
      FTLE_FILE_vtkDisplay.Origin = [0.0, 0.0, 0.0]
      FTLE_FILE_vtkDisplay.Pickable = 1
      FTLE_FILE_vtkDisplay.Texture = None
      FTLE_FILE_vtkDisplay.Triangulate = 0
      FTLE_FILE_vtkDisplay.UseShaderReplacements = 0
      FTLE_FILE_vtkDisplay.ShaderReplacements = ''
      FTLE_FILE_vtkDisplay.NonlinearSubdivisionLevel = 1
      FTLE_FILE_vtkDisplay.UseDataPartitions = 0
      FTLE_FILE_vtkDisplay.OSPRayUseScaleArray = 0
      FTLE_FILE_vtkDisplay.OSPRayScaleArray = 'FTLE - No T'
      FTLE_FILE_vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      FTLE_FILE_vtkDisplay.OSPRayMaterial = 'None'
      FTLE_FILE_vtkDisplay.Orient = 0
      FTLE_FILE_vtkDisplay.OrientationMode = 'Direction'
      FTLE_FILE_vtkDisplay.SelectOrientationVectors = 'FTLE - No T'
      FTLE_FILE_vtkDisplay.Scaling = 0
      FTLE_FILE_vtkDisplay.ScaleMode = 'No Data Scaling Off'
      FTLE_FILE_vtkDisplay.ScaleFactor = 4.5
      FTLE_FILE_vtkDisplay.SelectScaleArray = 'FTLE - No T'
      FTLE_FILE_vtkDisplay.GlyphType = 'Arrow'
      FTLE_FILE_vtkDisplay.UseGlyphTable = 0
      FTLE_FILE_vtkDisplay.GlyphTableIndexArray = 'FTLE - No T'
      FTLE_FILE_vtkDisplay.UseCompositeGlyphTable = 0
      FTLE_FILE_vtkDisplay.UseGlyphCullingAndLOD = 0
      FTLE_FILE_vtkDisplay.LODValues = []
      FTLE_FILE_vtkDisplay.ColorByLODIndex = 0
      FTLE_FILE_vtkDisplay.GaussianRadius = 0.225
      FTLE_FILE_vtkDisplay.ShaderPreset = 'Sphere'
      FTLE_FILE_vtkDisplay.CustomTriangleScale = 3
      FTLE_FILE_vtkDisplay.CustomShader = """ // This custom shader code define a gaussian blur
       // Please take a look into vtkSMPointGaussianRepresentation.cxx
       // for other custom shader examples
       //VTK::Color::Impl
         float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
         float gaussian = exp(-0.5*dist2);
         opacity = opacity*gaussian;
      """
      FTLE_FILE_vtkDisplay.Emissive = 0
      FTLE_FILE_vtkDisplay.ScaleByArray = 0
      FTLE_FILE_vtkDisplay.SetScaleArray = ['POINTS', 'FTLE - No T']
      FTLE_FILE_vtkDisplay.ScaleArrayComponent = ''
      FTLE_FILE_vtkDisplay.UseScaleFunction = 1
      FTLE_FILE_vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      FTLE_FILE_vtkDisplay.OpacityByArray = 0
      FTLE_FILE_vtkDisplay.OpacityArray = ['POINTS', 'FTLE - No T']
      FTLE_FILE_vtkDisplay.OpacityArrayComponent = ''
      FTLE_FILE_vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
      FTLE_FILE_vtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
      FTLE_FILE_vtkDisplay.SelectionCellLabelBold = 0
      FTLE_FILE_vtkDisplay.SelectionCellLabelColor = [0.0, 1.0, 0.0]
      FTLE_FILE_vtkDisplay.SelectionCellLabelFontFamily = 'Arial'
      FTLE_FILE_vtkDisplay.SelectionCellLabelFontFile = ''
      FTLE_FILE_vtkDisplay.SelectionCellLabelFontSize = 18
      FTLE_FILE_vtkDisplay.SelectionCellLabelItalic = 0
      FTLE_FILE_vtkDisplay.SelectionCellLabelJustification = 'Left'
      FTLE_FILE_vtkDisplay.SelectionCellLabelOpacity = 1.0
      FTLE_FILE_vtkDisplay.SelectionCellLabelShadow = 0
      FTLE_FILE_vtkDisplay.SelectionPointLabelBold = 0
      FTLE_FILE_vtkDisplay.SelectionPointLabelColor = [1.0, 1.0, 0.0]
      FTLE_FILE_vtkDisplay.SelectionPointLabelFontFamily = 'Arial'
      FTLE_FILE_vtkDisplay.SelectionPointLabelFontFile = ''
      FTLE_FILE_vtkDisplay.SelectionPointLabelFontSize = 18
      FTLE_FILE_vtkDisplay.SelectionPointLabelItalic = 0
      FTLE_FILE_vtkDisplay.SelectionPointLabelJustification = 'Left'
      FTLE_FILE_vtkDisplay.SelectionPointLabelOpacity = 1.0
      FTLE_FILE_vtkDisplay.SelectionPointLabelShadow = 0
      FTLE_FILE_vtkDisplay.PolarAxes = 'PolarAxesRepresentation'
      FTLE_FILE_vtkDisplay.ScalarOpacityFunction = None
      FTLE_FILE_vtkDisplay.ScalarOpacityUnitDistance = 0.5175326316478011
      FTLE_FILE_vtkDisplay.SelectMapper = 'Projected tetra'
      FTLE_FILE_vtkDisplay.SamplingDimensions = [128, 128, 128]

      # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
      FTLE_FILE_vtkDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
      FTLE_FILE_vtkDisplay.OSPRayScaleFunction.UseLogScale = 0

      # init the 'Arrow' selected for 'GlyphType'
      FTLE_FILE_vtkDisplay.GlyphType.TipResolution = 6
      FTLE_FILE_vtkDisplay.GlyphType.TipRadius = 0.1
      FTLE_FILE_vtkDisplay.GlyphType.TipLength = 0.35
      FTLE_FILE_vtkDisplay.GlyphType.ShaftResolution = 6
      FTLE_FILE_vtkDisplay.GlyphType.ShaftRadius = 0.03
      FTLE_FILE_vtkDisplay.GlyphType.Invert = 0

      # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
      FTLE_FILE_vtkDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 6.1694401555, 1.0, 0.5, 0.0]
      FTLE_FILE_vtkDisplay.ScaleTransferFunction.UseLogScale = 0

      # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
      FTLE_FILE_vtkDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 6.1694401555, 1.0, 0.5, 0.0]
      FTLE_FILE_vtkDisplay.OpacityTransferFunction.UseLogScale = 0

      # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
      FTLE_FILE_vtkDisplay.DataAxesGrid.XTitle = 'X Axis'
      FTLE_FILE_vtkDisplay.DataAxesGrid.YTitle = 'Y Axis'
      FTLE_FILE_vtkDisplay.DataAxesGrid.ZTitle = 'Z Axis'
      FTLE_FILE_vtkDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
      FTLE_FILE_vtkDisplay.DataAxesGrid.XTitleFontFamily = 'Arial'
      FTLE_FILE_vtkDisplay.DataAxesGrid.XTitleFontFile = ''
      FTLE_FILE_vtkDisplay.DataAxesGrid.XTitleBold = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.XTitleItalic = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.XTitleFontSize = 12
      FTLE_FILE_vtkDisplay.DataAxesGrid.XTitleShadow = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.XTitleOpacity = 1.0
      FTLE_FILE_vtkDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
      FTLE_FILE_vtkDisplay.DataAxesGrid.YTitleFontFamily = 'Arial'
      FTLE_FILE_vtkDisplay.DataAxesGrid.YTitleFontFile = ''
      FTLE_FILE_vtkDisplay.DataAxesGrid.YTitleBold = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.YTitleItalic = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.YTitleFontSize = 12
      FTLE_FILE_vtkDisplay.DataAxesGrid.YTitleShadow = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.YTitleOpacity = 1.0
      FTLE_FILE_vtkDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
      FTLE_FILE_vtkDisplay.DataAxesGrid.ZTitleFontFamily = 'Arial'
      FTLE_FILE_vtkDisplay.DataAxesGrid.ZTitleFontFile = ''
      FTLE_FILE_vtkDisplay.DataAxesGrid.ZTitleBold = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.ZTitleItalic = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.ZTitleFontSize = 12
      FTLE_FILE_vtkDisplay.DataAxesGrid.ZTitleShadow = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.ZTitleOpacity = 1.0
      FTLE_FILE_vtkDisplay.DataAxesGrid.FacesToRender = 63
      FTLE_FILE_vtkDisplay.DataAxesGrid.CullBackface = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.CullFrontface = 1
      FTLE_FILE_vtkDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
      FTLE_FILE_vtkDisplay.DataAxesGrid.ShowGrid = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.ShowEdges = 1
      FTLE_FILE_vtkDisplay.DataAxesGrid.ShowTicks = 1
      FTLE_FILE_vtkDisplay.DataAxesGrid.LabelUniqueEdgesOnly = 1
      FTLE_FILE_vtkDisplay.DataAxesGrid.AxesToLabel = 63
      FTLE_FILE_vtkDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
      FTLE_FILE_vtkDisplay.DataAxesGrid.XLabelFontFamily = 'Arial'
      FTLE_FILE_vtkDisplay.DataAxesGrid.XLabelFontFile = ''
      FTLE_FILE_vtkDisplay.DataAxesGrid.XLabelBold = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.XLabelItalic = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.XLabelFontSize = 12
      FTLE_FILE_vtkDisplay.DataAxesGrid.XLabelShadow = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.XLabelOpacity = 1.0
      FTLE_FILE_vtkDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
      FTLE_FILE_vtkDisplay.DataAxesGrid.YLabelFontFamily = 'Arial'
      FTLE_FILE_vtkDisplay.DataAxesGrid.YLabelFontFile = ''
      FTLE_FILE_vtkDisplay.DataAxesGrid.YLabelBold = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.YLabelItalic = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.YLabelFontSize = 12
      FTLE_FILE_vtkDisplay.DataAxesGrid.YLabelShadow = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.YLabelOpacity = 1.0
      FTLE_FILE_vtkDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
      FTLE_FILE_vtkDisplay.DataAxesGrid.ZLabelFontFamily = 'Arial'
      FTLE_FILE_vtkDisplay.DataAxesGrid.ZLabelFontFile = ''
      FTLE_FILE_vtkDisplay.DataAxesGrid.ZLabelBold = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.ZLabelItalic = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.ZLabelFontSize = 12
      FTLE_FILE_vtkDisplay.DataAxesGrid.ZLabelShadow = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.ZLabelOpacity = 1.0
      FTLE_FILE_vtkDisplay.DataAxesGrid.XAxisNotation = 'Mixed'
      FTLE_FILE_vtkDisplay.DataAxesGrid.XAxisPrecision = 2
      FTLE_FILE_vtkDisplay.DataAxesGrid.XAxisUseCustomLabels = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.XAxisLabels = []
      FTLE_FILE_vtkDisplay.DataAxesGrid.YAxisNotation = 'Mixed'
      FTLE_FILE_vtkDisplay.DataAxesGrid.YAxisPrecision = 2
      FTLE_FILE_vtkDisplay.DataAxesGrid.YAxisUseCustomLabels = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.YAxisLabels = []
      FTLE_FILE_vtkDisplay.DataAxesGrid.ZAxisNotation = 'Mixed'
      FTLE_FILE_vtkDisplay.DataAxesGrid.ZAxisPrecision = 2
      FTLE_FILE_vtkDisplay.DataAxesGrid.ZAxisUseCustomLabels = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.ZAxisLabels = []
      FTLE_FILE_vtkDisplay.DataAxesGrid.UseCustomBounds = 0
      FTLE_FILE_vtkDisplay.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

      # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
      FTLE_FILE_vtkDisplay.PolarAxes.Visibility = 0
      FTLE_FILE_vtkDisplay.PolarAxes.Translation = [0.0, 0.0, 0.0]
      FTLE_FILE_vtkDisplay.PolarAxes.Scale = [1.0, 1.0, 1.0]
      FTLE_FILE_vtkDisplay.PolarAxes.Orientation = [0.0, 0.0, 0.0]
      FTLE_FILE_vtkDisplay.PolarAxes.EnableCustomBounds = [0, 0, 0]
      FTLE_FILE_vtkDisplay.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
      FTLE_FILE_vtkDisplay.PolarAxes.EnableCustomRange = 0
      FTLE_FILE_vtkDisplay.PolarAxes.CustomRange = [0.0, 1.0]
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisVisibility = 1
      FTLE_FILE_vtkDisplay.PolarAxes.RadialAxesVisibility = 1
      FTLE_FILE_vtkDisplay.PolarAxes.DrawRadialGridlines = 1
      FTLE_FILE_vtkDisplay.PolarAxes.PolarArcsVisibility = 1
      FTLE_FILE_vtkDisplay.PolarAxes.DrawPolarArcsGridlines = 1
      FTLE_FILE_vtkDisplay.PolarAxes.NumberOfRadialAxes = 0
      FTLE_FILE_vtkDisplay.PolarAxes.AutoSubdividePolarAxis = 1
      FTLE_FILE_vtkDisplay.PolarAxes.NumberOfPolarAxis = 0
      FTLE_FILE_vtkDisplay.PolarAxes.MinimumRadius = 0.0
      FTLE_FILE_vtkDisplay.PolarAxes.MinimumAngle = 0.0
      FTLE_FILE_vtkDisplay.PolarAxes.MaximumAngle = 90.0
      FTLE_FILE_vtkDisplay.PolarAxes.RadialAxesOriginToPolarAxis = 1
      FTLE_FILE_vtkDisplay.PolarAxes.Ratio = 1.0
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
      FTLE_FILE_vtkDisplay.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
      FTLE_FILE_vtkDisplay.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
      FTLE_FILE_vtkDisplay.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
      FTLE_FILE_vtkDisplay.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisTitleVisibility = 1
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisTitle = 'Radial Distance'
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisTitleLocation = 'Bottom'
      FTLE_FILE_vtkDisplay.PolarAxes.PolarLabelVisibility = 1
      FTLE_FILE_vtkDisplay.PolarAxes.PolarLabelFormat = '%-#6.3g'
      FTLE_FILE_vtkDisplay.PolarAxes.PolarLabelExponentLocation = 'Labels'
      FTLE_FILE_vtkDisplay.PolarAxes.RadialLabelVisibility = 1
      FTLE_FILE_vtkDisplay.PolarAxes.RadialLabelFormat = '%-#3.1f'
      FTLE_FILE_vtkDisplay.PolarAxes.RadialLabelLocation = 'Bottom'
      FTLE_FILE_vtkDisplay.PolarAxes.RadialUnitsVisibility = 1
      FTLE_FILE_vtkDisplay.PolarAxes.ScreenSize = 10.0
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisTitleOpacity = 1.0
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisTitleFontFile = ''
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisTitleBold = 0
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisTitleItalic = 0
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisTitleShadow = 0
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisTitleFontSize = 12
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisLabelOpacity = 1.0
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisLabelFontFile = ''
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisLabelBold = 0
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisLabelItalic = 0
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisLabelShadow = 0
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisLabelFontSize = 12
      FTLE_FILE_vtkDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
      FTLE_FILE_vtkDisplay.PolarAxes.LastRadialAxisTextOpacity = 1.0
      FTLE_FILE_vtkDisplay.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
      FTLE_FILE_vtkDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
      FTLE_FILE_vtkDisplay.PolarAxes.LastRadialAxisTextBold = 0
      FTLE_FILE_vtkDisplay.PolarAxes.LastRadialAxisTextItalic = 0
      FTLE_FILE_vtkDisplay.PolarAxes.LastRadialAxisTextShadow = 0
      FTLE_FILE_vtkDisplay.PolarAxes.LastRadialAxisTextFontSize = 12
      FTLE_FILE_vtkDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
      FTLE_FILE_vtkDisplay.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
      FTLE_FILE_vtkDisplay.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
      FTLE_FILE_vtkDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''
      FTLE_FILE_vtkDisplay.PolarAxes.SecondaryRadialAxesTextBold = 0
      FTLE_FILE_vtkDisplay.PolarAxes.SecondaryRadialAxesTextItalic = 0
      FTLE_FILE_vtkDisplay.PolarAxes.SecondaryRadialAxesTextShadow = 0
      FTLE_FILE_vtkDisplay.PolarAxes.SecondaryRadialAxesTextFontSize = 12
      FTLE_FILE_vtkDisplay.PolarAxes.EnableDistanceLOD = 1
      FTLE_FILE_vtkDisplay.PolarAxes.DistanceLODThreshold = 0.7
      FTLE_FILE_vtkDisplay.PolarAxes.EnableViewAngleLOD = 1
      FTLE_FILE_vtkDisplay.PolarAxes.ViewAngleLODThreshold = 0.7
      FTLE_FILE_vtkDisplay.PolarAxes.SmallestVisiblePolarAngle = 0.5
      FTLE_FILE_vtkDisplay.PolarAxes.PolarTicksVisibility = 1
      FTLE_FILE_vtkDisplay.PolarAxes.ArcTicksOriginToPolarAxis = 1
      FTLE_FILE_vtkDisplay.PolarAxes.TickLocation = 'Both'
      FTLE_FILE_vtkDisplay.PolarAxes.AxisTickVisibility = 1
      FTLE_FILE_vtkDisplay.PolarAxes.AxisMinorTickVisibility = 0
      FTLE_FILE_vtkDisplay.PolarAxes.ArcTickVisibility = 1
      FTLE_FILE_vtkDisplay.PolarAxes.ArcMinorTickVisibility = 0
      FTLE_FILE_vtkDisplay.PolarAxes.DeltaAngleMajor = 10.0
      FTLE_FILE_vtkDisplay.PolarAxes.DeltaAngleMinor = 5.0
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisMajorTickSize = 0.0
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisTickRatioSize = 0.3
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisMajorTickThickness = 1.0
      FTLE_FILE_vtkDisplay.PolarAxes.PolarAxisTickRatioThickness = 0.5
      FTLE_FILE_vtkDisplay.PolarAxes.LastRadialAxisMajorTickSize = 0.0
      FTLE_FILE_vtkDisplay.PolarAxes.LastRadialAxisTickRatioSize = 0.3
      FTLE_FILE_vtkDisplay.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
      FTLE_FILE_vtkDisplay.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
      FTLE_FILE_vtkDisplay.PolarAxes.ArcMajorTickSize = 0.0
      FTLE_FILE_vtkDisplay.PolarAxes.ArcTickRatioSize = 0.3
      FTLE_FILE_vtkDisplay.PolarAxes.ArcMajorTickThickness = 1.0
      FTLE_FILE_vtkDisplay.PolarAxes.ArcTickRatioThickness = 0.5
      FTLE_FILE_vtkDisplay.PolarAxes.Use2DMode = 0
      FTLE_FILE_vtkDisplay.PolarAxes.UseLogAxis = 0

      # reset view to fit data
      renderView1.ResetCamera()

      #changing interaction mode based on data extents
      renderView1.InteractionMode = '2D'
      renderView1.CameraPosition = [22.5, 3.0, 10000.0]
      renderView1.CameraFocalPoint = [22.5, 3.0, 0.0]

      # update the view to ensure updated data information
      renderView1.Update()

      # set scalar coloring
      ColorBy(FTLE_FILE_vtkDisplay, ('POINTS', 'FTLE - No T'))

      # rescale color and/or opacity maps used to include current data range
      FTLE_FILE_vtkDisplay.RescaleTransferFunctionToDataRange(True, False)

      # show color bar/color legend
      FTLE_FILE_vtkDisplay.SetScalarBarVisibility(renderView1, True)

      # get color transfer function/color map for 'FTLENoT'
      fTLENoTLUT = GetColorTransferFunction('FTLENoT')
      fTLENoTLUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
      fTLENoTLUT.InterpretValuesAsCategories = 0
      fTLENoTLUT.AnnotationsInitialized = 0
      fTLENoTLUT.ShowCategoricalColorsinDataRangeOnly = 0
      fTLENoTLUT.RescaleOnVisibilityChange = 0
      fTLENoTLUT.EnableOpacityMapping = 0
      fTLENoTLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 3.08472007775, 0.865003, 0.865003, 0.865003, 6.1694401555, 0.705882, 0.0156863, 0.14902]
      fTLENoTLUT.UseLogScale = 0
      fTLENoTLUT.ColorSpace = 'Diverging'
      fTLENoTLUT.UseBelowRangeColor = 0
      fTLENoTLUT.BelowRangeColor = [0.0, 0.0, 0.0]
      fTLENoTLUT.UseAboveRangeColor = 0
      fTLENoTLUT.AboveRangeColor = [0.5, 0.5, 0.5]
      fTLENoTLUT.NanColor = [1.0, 1.0, 0.0]
      fTLENoTLUT.NanOpacity = 1.0
      fTLENoTLUT.Discretize = 1
      fTLENoTLUT.NumberOfTableValues = 256
      fTLENoTLUT.ScalarRangeInitialized = 1.0
      fTLENoTLUT.HSVWrap = 0
      fTLENoTLUT.VectorComponent = 0
      fTLENoTLUT.VectorMode = 'Magnitude'
      fTLENoTLUT.AllowDuplicateScalars = 1
      fTLENoTLUT.Annotations = []
      fTLENoTLUT.ActiveAnnotatedValues = []
      fTLENoTLUT.IndexedColors = []
      fTLENoTLUT.IndexedOpacities = []

      # get opacity transfer function/opacity map for 'FTLENoT'
      fTLENoTPWF = GetOpacityTransferFunction('FTLENoT')
      fTLENoTPWF.Points = [0.0, 0.0, 0.5, 0.0, 6.1694401555, 1.0, 0.5, 0.0]
      fTLENoTPWF.AllowDuplicateScalars = 1
      fTLENoTPWF.UseLogScale = 0
      fTLENoTPWF.ScalarRangeInitialized = 1

      # hide color bar/color legend
      FTLE_FILE_vtkDisplay.SetScalarBarVisibility(renderView1, False)

      # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
      fTLENoTLUT.ApplyPreset('jet', True)

      # Rescale transfer function
      fTLENoTLUT.RescaleTransferFunction(0.0, 5.0)

      # Rescale transfer function
      fTLENoTPWF.RescaleTransferFunction(0.0, 5.0)

      # update the view to ensure updated data information
      renderView1.Update()

      # current camera placement for renderView1
      renderView1.InteractionMode = '2D'
      renderView1.CameraPosition = [22.5, 3.0, 10000.0]
      renderView1.CameraFocalPoint = [22.5, 3.0, 0.0]
      renderView1.CameraParallelScale = 14

      # save screenshot
      SaveScreenshot(outName, renderView1, ImageResolution=[4096, 2160],
          FontScaling='Scale fonts proportionally',
          OverrideColorPalette='',
          StereoMode='No change',
          TransparentBackground=0, 
          # PNG options
          CompressionLevel='5')

      # set active source
      SetActiveSource(FTLE_FILE_vtk)

      # destroy FTLE_FILE_vtk
      Delete(FTLE_FILE_vtk)
      del FTLE_FILE_vtk

      #### saving camera placements for all active views

      # current camera placement for renderView1
      renderView1.InteractionMode = '2D'
      renderView1.CameraPosition = [22.5, 3.0, 10000.0]
      renderView1.CameraFocalPoint = [22.5, 3.0, 0.0]
      renderView1.CameraParallelScale = 12.813060868709753

      print('Wrote screenshot to %s' %outName)
