import sys, os
from paraview.simple import *

structNames = ['M1', 'M2', 'M3', 'M4', 'M7', 'M8', 'M9']
leakNames   = ['0PL', '10PL', '20PL', '40PL']
timeIndexes = ['6','11','14','24','30','38','50','68']
timeNames   = ['T1','T2','T3','T4','T5','T6','T7','T8']

baseClotDir = os.getcwd() + '/Clot-Microstructures/'
baseFTLEDir = os.getcwd() + '/output/'
FTLEDir = '/media/zach/BACKUP/Mukherjee/output/Paper/LS - T=0.05/'

for struct in structNames:

  structFileName = baseClotDir + 'clot-geo-high-shear-final.vtk'

  for leak in leakNames:

    leakNo = leak.split('P')[0]

    for timeInd in range(len(timeIndexes)):

      ftleFileName = baseFTLEDir + 'FTLE-HS-%s-%s/FTLE-HS-%s-%spL-%s.vtk' %(struct, leak, struct, leakNo, timeIndexes[timeInd])

      outName      = baseFTLEDir + 'FTLE_Geo/High-Shear-Clot/%s/%s-Percent-Leakage/FTLE-HS-%s-%s-%s.png' %(struct, leakNo, struct, leakNo, timeNames[timeInd])

      paraview.simple._DisableFirstRenderCameraReset()

      LoadPalette(paletteName='WhiteBackground')

      # get active view
      renderView1 = GetActiveViewOrCreate('RenderView')

      renderView1.UseLight = 0

      # Properties modified on renderView1
      renderView1.OrientationAxesVisibility = 0

      # get the material library
      materialLibrary1 = GetMaterialLibrary()

      # create a new 'Legacy VTK Reader'
      FTLE_FILE = LegacyVTKReader(FileNames=[ftleFileName])

      # show data in view
      FTLE_FILEDisplay = Show(FTLE_FILE, renderView1)

      # trace defaults for the display properties.
      FTLE_FILEDisplay.Representation = 'Surface'
      FTLE_FILEDisplay.AmbientColor = [0.0, 0.0, 0.0]
      FTLE_FILEDisplay.ColorArrayName = [None, '']
      FTLE_FILEDisplay.DiffuseColor = [1.0, 1.0, 1.0]
      FTLE_FILEDisplay.LookupTable = None
      FTLE_FILEDisplay.MapScalars = 1
      FTLE_FILEDisplay.MultiComponentsMapping = 0
      FTLE_FILEDisplay.InterpolateScalarsBeforeMapping = 1
      FTLE_FILEDisplay.Opacity = 1.0
      FTLE_FILEDisplay.PointSize = 2.0
      FTLE_FILEDisplay.LineWidth = 1.0
      FTLE_FILEDisplay.RenderLinesAsTubes = 0
      FTLE_FILEDisplay.RenderPointsAsSpheres = 0
      FTLE_FILEDisplay.Interpolation = 'Gouraud'
      FTLE_FILEDisplay.Specular = 0.0
      FTLE_FILEDisplay.SpecularColor = [1.0, 1.0, 1.0]
      FTLE_FILEDisplay.SpecularPower = 100.0
      FTLE_FILEDisplay.Luminosity = 0.0
      FTLE_FILEDisplay.Ambient = 0.0
      FTLE_FILEDisplay.Diffuse = 1.0
      FTLE_FILEDisplay.EdgeColor = [0.0, 0.0, 0.5]
      FTLE_FILEDisplay.BackfaceRepresentation = 'Follow Frontface'
      FTLE_FILEDisplay.BackfaceAmbientColor = [1.0, 1.0, 1.0]
      FTLE_FILEDisplay.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
      FTLE_FILEDisplay.BackfaceOpacity = 1.0
      FTLE_FILEDisplay.Position = [0.0, 0.0, 0.0]
      FTLE_FILEDisplay.Scale = [1.0, 1.0, 1.0]
      FTLE_FILEDisplay.Orientation = [0.0, 0.0, 0.0]
      FTLE_FILEDisplay.Origin = [0.0, 0.0, 0.0]
      FTLE_FILEDisplay.Pickable = 1
      FTLE_FILEDisplay.Texture = None
      FTLE_FILEDisplay.Triangulate = 0
      FTLE_FILEDisplay.UseShaderReplacements = 0
      FTLE_FILEDisplay.ShaderReplacements = ''
      FTLE_FILEDisplay.NonlinearSubdivisionLevel = 1
      FTLE_FILEDisplay.UseDataPartitions = 0
      FTLE_FILEDisplay.OSPRayUseScaleArray = 0
      FTLE_FILEDisplay.OSPRayScaleArray = 'FTLE - No T'
      FTLE_FILEDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      FTLE_FILEDisplay.OSPRayMaterial = 'None'
      FTLE_FILEDisplay.Orient = 0
      FTLE_FILEDisplay.OrientationMode = 'Direction'
      FTLE_FILEDisplay.SelectOrientationVectors = 'FTLE - No T'
      FTLE_FILEDisplay.Scaling = 0
      FTLE_FILEDisplay.ScaleMode = 'No Data Scaling Off'
      FTLE_FILEDisplay.ScaleFactor = 4.5
      FTLE_FILEDisplay.SelectScaleArray = 'FTLE - No T'
      FTLE_FILEDisplay.GlyphType = 'Arrow'
      FTLE_FILEDisplay.UseGlyphTable = 0
      FTLE_FILEDisplay.GlyphTableIndexArray = 'FTLE - No T'
      FTLE_FILEDisplay.UseCompositeGlyphTable = 0
      FTLE_FILEDisplay.UseGlyphCullingAndLOD = 0
      FTLE_FILEDisplay.LODValues = []
      FTLE_FILEDisplay.ColorByLODIndex = 0
      FTLE_FILEDisplay.GaussianRadius = 0.225
      FTLE_FILEDisplay.ShaderPreset = 'Sphere'
      FTLE_FILEDisplay.CustomTriangleScale = 3
      FTLE_FILEDisplay.CustomShader = """ // This custom shader code define a gaussian blur
       // Please take a look into vtkSMPointGaussianRepresentation.cxx
       // for other custom shader examples
       //VTK::Color::Impl
         float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
         float gaussian = exp(-0.5*dist2);
         opacity = opacity*gaussian;
      """
      FTLE_FILEDisplay.Emissive = 0
      FTLE_FILEDisplay.ScaleByArray = 0
      FTLE_FILEDisplay.SetScaleArray = ['POINTS', 'FTLE - No T']
      FTLE_FILEDisplay.ScaleArrayComponent = ''
      FTLE_FILEDisplay.UseScaleFunction = 1
      FTLE_FILEDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      FTLE_FILEDisplay.OpacityByArray = 0
      FTLE_FILEDisplay.OpacityArray = ['POINTS', 'FTLE - No T']
      FTLE_FILEDisplay.OpacityArrayComponent = ''
      FTLE_FILEDisplay.OpacityTransferFunction = 'PiecewiseFunction'
      FTLE_FILEDisplay.DataAxesGrid = 'GridAxesRepresentation'
      FTLE_FILEDisplay.SelectionCellLabelBold = 0
      FTLE_FILEDisplay.SelectionCellLabelColor = [0.0, 1.0, 0.0]
      FTLE_FILEDisplay.SelectionCellLabelFontFamily = 'Arial'
      FTLE_FILEDisplay.SelectionCellLabelFontFile = ''
      FTLE_FILEDisplay.SelectionCellLabelFontSize = 18
      FTLE_FILEDisplay.SelectionCellLabelItalic = 0
      FTLE_FILEDisplay.SelectionCellLabelJustification = 'Left'
      FTLE_FILEDisplay.SelectionCellLabelOpacity = 1.0
      FTLE_FILEDisplay.SelectionCellLabelShadow = 0
      FTLE_FILEDisplay.SelectionPointLabelBold = 0
      FTLE_FILEDisplay.SelectionPointLabelColor = [1.0, 1.0, 0.0]
      FTLE_FILEDisplay.SelectionPointLabelFontFamily = 'Arial'
      FTLE_FILEDisplay.SelectionPointLabelFontFile = ''
      FTLE_FILEDisplay.SelectionPointLabelFontSize = 18
      FTLE_FILEDisplay.SelectionPointLabelItalic = 0
      FTLE_FILEDisplay.SelectionPointLabelJustification = 'Left'
      FTLE_FILEDisplay.SelectionPointLabelOpacity = 1.0
      FTLE_FILEDisplay.SelectionPointLabelShadow = 0
      FTLE_FILEDisplay.PolarAxes = 'PolarAxesRepresentation'
      FTLE_FILEDisplay.ScalarOpacityFunction = None
      FTLE_FILEDisplay.ScalarOpacityUnitDistance = 0.5175326316478011
      FTLE_FILEDisplay.SelectMapper = 'Projected tetra'
      FTLE_FILEDisplay.SamplingDimensions = [128, 128, 128]

      # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
      FTLE_FILEDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
      FTLE_FILEDisplay.OSPRayScaleFunction.UseLogScale = 0

      # init the 'Arrow' selected for 'GlyphType'
      FTLE_FILEDisplay.GlyphType.TipResolution = 6
      FTLE_FILEDisplay.GlyphType.TipRadius = 0.1
      FTLE_FILEDisplay.GlyphType.TipLength = 0.35
      FTLE_FILEDisplay.GlyphType.ShaftResolution = 6
      FTLE_FILEDisplay.GlyphType.ShaftRadius = 0.03
      FTLE_FILEDisplay.GlyphType.Invert = 0

      # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
      FTLE_FILEDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 4.879917279551101, 1.0, 0.5, 0.0]
      FTLE_FILEDisplay.ScaleTransferFunction.UseLogScale = 0

      # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
      FTLE_FILEDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 4.879917279551101, 1.0, 0.5, 0.0]
      FTLE_FILEDisplay.OpacityTransferFunction.UseLogScale = 0

      # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
      FTLE_FILEDisplay.DataAxesGrid.XTitle = 'X Axis'
      FTLE_FILEDisplay.DataAxesGrid.YTitle = 'Y Axis'
      FTLE_FILEDisplay.DataAxesGrid.ZTitle = 'Z Axis'
      FTLE_FILEDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
      FTLE_FILEDisplay.DataAxesGrid.XTitleFontFamily = 'Arial'
      FTLE_FILEDisplay.DataAxesGrid.XTitleFontFile = ''
      FTLE_FILEDisplay.DataAxesGrid.XTitleBold = 0
      FTLE_FILEDisplay.DataAxesGrid.XTitleItalic = 0
      FTLE_FILEDisplay.DataAxesGrid.XTitleFontSize = 12
      FTLE_FILEDisplay.DataAxesGrid.XTitleShadow = 0
      FTLE_FILEDisplay.DataAxesGrid.XTitleOpacity = 1.0
      FTLE_FILEDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
      FTLE_FILEDisplay.DataAxesGrid.YTitleFontFamily = 'Arial'
      FTLE_FILEDisplay.DataAxesGrid.YTitleFontFile = ''
      FTLE_FILEDisplay.DataAxesGrid.YTitleBold = 0
      FTLE_FILEDisplay.DataAxesGrid.YTitleItalic = 0
      FTLE_FILEDisplay.DataAxesGrid.YTitleFontSize = 12
      FTLE_FILEDisplay.DataAxesGrid.YTitleShadow = 0
      FTLE_FILEDisplay.DataAxesGrid.YTitleOpacity = 1.0
      FTLE_FILEDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
      FTLE_FILEDisplay.DataAxesGrid.ZTitleFontFamily = 'Arial'
      FTLE_FILEDisplay.DataAxesGrid.ZTitleFontFile = ''
      FTLE_FILEDisplay.DataAxesGrid.ZTitleBold = 0
      FTLE_FILEDisplay.DataAxesGrid.ZTitleItalic = 0
      FTLE_FILEDisplay.DataAxesGrid.ZTitleFontSize = 12
      FTLE_FILEDisplay.DataAxesGrid.ZTitleShadow = 0
      FTLE_FILEDisplay.DataAxesGrid.ZTitleOpacity = 1.0
      FTLE_FILEDisplay.DataAxesGrid.FacesToRender = 63
      FTLE_FILEDisplay.DataAxesGrid.CullBackface = 0
      FTLE_FILEDisplay.DataAxesGrid.CullFrontface = 1
      FTLE_FILEDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
      FTLE_FILEDisplay.DataAxesGrid.ShowGrid = 0
      FTLE_FILEDisplay.DataAxesGrid.ShowEdges = 1
      FTLE_FILEDisplay.DataAxesGrid.ShowTicks = 1
      FTLE_FILEDisplay.DataAxesGrid.LabelUniqueEdgesOnly = 1
      FTLE_FILEDisplay.DataAxesGrid.AxesToLabel = 63
      FTLE_FILEDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
      FTLE_FILEDisplay.DataAxesGrid.XLabelFontFamily = 'Arial'
      FTLE_FILEDisplay.DataAxesGrid.XLabelFontFile = ''
      FTLE_FILEDisplay.DataAxesGrid.XLabelBold = 0
      FTLE_FILEDisplay.DataAxesGrid.XLabelItalic = 0
      FTLE_FILEDisplay.DataAxesGrid.XLabelFontSize = 12
      FTLE_FILEDisplay.DataAxesGrid.XLabelShadow = 0
      FTLE_FILEDisplay.DataAxesGrid.XLabelOpacity = 1.0
      FTLE_FILEDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
      FTLE_FILEDisplay.DataAxesGrid.YLabelFontFamily = 'Arial'
      FTLE_FILEDisplay.DataAxesGrid.YLabelFontFile = ''
      FTLE_FILEDisplay.DataAxesGrid.YLabelBold = 0
      FTLE_FILEDisplay.DataAxesGrid.YLabelItalic = 0
      FTLE_FILEDisplay.DataAxesGrid.YLabelFontSize = 12
      FTLE_FILEDisplay.DataAxesGrid.YLabelShadow = 0
      FTLE_FILEDisplay.DataAxesGrid.YLabelOpacity = 1.0
      FTLE_FILEDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]
      FTLE_FILEDisplay.DataAxesGrid.ZLabelFontFamily = 'Arial'
      FTLE_FILEDisplay.DataAxesGrid.ZLabelFontFile = ''
      FTLE_FILEDisplay.DataAxesGrid.ZLabelBold = 0
      FTLE_FILEDisplay.DataAxesGrid.ZLabelItalic = 0
      FTLE_FILEDisplay.DataAxesGrid.ZLabelFontSize = 12
      FTLE_FILEDisplay.DataAxesGrid.ZLabelShadow = 0
      FTLE_FILEDisplay.DataAxesGrid.ZLabelOpacity = 1.0
      FTLE_FILEDisplay.DataAxesGrid.XAxisNotation = 'Mixed'
      FTLE_FILEDisplay.DataAxesGrid.XAxisPrecision = 2
      FTLE_FILEDisplay.DataAxesGrid.XAxisUseCustomLabels = 0
      FTLE_FILEDisplay.DataAxesGrid.XAxisLabels = []
      FTLE_FILEDisplay.DataAxesGrid.YAxisNotation = 'Mixed'
      FTLE_FILEDisplay.DataAxesGrid.YAxisPrecision = 2
      FTLE_FILEDisplay.DataAxesGrid.YAxisUseCustomLabels = 0
      FTLE_FILEDisplay.DataAxesGrid.YAxisLabels = []
      FTLE_FILEDisplay.DataAxesGrid.ZAxisNotation = 'Mixed'
      FTLE_FILEDisplay.DataAxesGrid.ZAxisPrecision = 2
      FTLE_FILEDisplay.DataAxesGrid.ZAxisUseCustomLabels = 0
      FTLE_FILEDisplay.DataAxesGrid.ZAxisLabels = []
      FTLE_FILEDisplay.DataAxesGrid.UseCustomBounds = 0
      FTLE_FILEDisplay.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

      # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
      FTLE_FILEDisplay.PolarAxes.Visibility = 0
      FTLE_FILEDisplay.PolarAxes.Translation = [0.0, 0.0, 0.0]
      FTLE_FILEDisplay.PolarAxes.Scale = [1.0, 1.0, 1.0]
      FTLE_FILEDisplay.PolarAxes.Orientation = [0.0, 0.0, 0.0]
      FTLE_FILEDisplay.PolarAxes.EnableCustomBounds = [0, 0, 0]
      FTLE_FILEDisplay.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
      FTLE_FILEDisplay.PolarAxes.EnableCustomRange = 0
      FTLE_FILEDisplay.PolarAxes.CustomRange = [0.0, 1.0]
      FTLE_FILEDisplay.PolarAxes.PolarAxisVisibility = 1
      FTLE_FILEDisplay.PolarAxes.RadialAxesVisibility = 1
      FTLE_FILEDisplay.PolarAxes.DrawRadialGridlines = 1
      FTLE_FILEDisplay.PolarAxes.PolarArcsVisibility = 1
      FTLE_FILEDisplay.PolarAxes.DrawPolarArcsGridlines = 1
      FTLE_FILEDisplay.PolarAxes.NumberOfRadialAxes = 0
      FTLE_FILEDisplay.PolarAxes.AutoSubdividePolarAxis = 1
      FTLE_FILEDisplay.PolarAxes.NumberOfPolarAxis = 0
      FTLE_FILEDisplay.PolarAxes.MinimumRadius = 0.0
      FTLE_FILEDisplay.PolarAxes.MinimumAngle = 0.0
      FTLE_FILEDisplay.PolarAxes.MaximumAngle = 90.0
      FTLE_FILEDisplay.PolarAxes.RadialAxesOriginToPolarAxis = 1
      FTLE_FILEDisplay.PolarAxes.Ratio = 1.0
      FTLE_FILEDisplay.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
      FTLE_FILEDisplay.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
      FTLE_FILEDisplay.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
      FTLE_FILEDisplay.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
      FTLE_FILEDisplay.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
      FTLE_FILEDisplay.PolarAxes.PolarAxisTitleVisibility = 1
      FTLE_FILEDisplay.PolarAxes.PolarAxisTitle = 'Radial Distance'
      FTLE_FILEDisplay.PolarAxes.PolarAxisTitleLocation = 'Bottom'
      FTLE_FILEDisplay.PolarAxes.PolarLabelVisibility = 1
      FTLE_FILEDisplay.PolarAxes.PolarLabelFormat = '%-#6.3g'
      FTLE_FILEDisplay.PolarAxes.PolarLabelExponentLocation = 'Labels'
      FTLE_FILEDisplay.PolarAxes.RadialLabelVisibility = 1
      FTLE_FILEDisplay.PolarAxes.RadialLabelFormat = '%-#3.1f'
      FTLE_FILEDisplay.PolarAxes.RadialLabelLocation = 'Bottom'
      FTLE_FILEDisplay.PolarAxes.RadialUnitsVisibility = 1
      FTLE_FILEDisplay.PolarAxes.ScreenSize = 10.0
      FTLE_FILEDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
      FTLE_FILEDisplay.PolarAxes.PolarAxisTitleOpacity = 1.0
      FTLE_FILEDisplay.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
      FTLE_FILEDisplay.PolarAxes.PolarAxisTitleFontFile = ''
      FTLE_FILEDisplay.PolarAxes.PolarAxisTitleBold = 0
      FTLE_FILEDisplay.PolarAxes.PolarAxisTitleItalic = 0
      FTLE_FILEDisplay.PolarAxes.PolarAxisTitleShadow = 0
      FTLE_FILEDisplay.PolarAxes.PolarAxisTitleFontSize = 12
      FTLE_FILEDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
      FTLE_FILEDisplay.PolarAxes.PolarAxisLabelOpacity = 1.0
      FTLE_FILEDisplay.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
      FTLE_FILEDisplay.PolarAxes.PolarAxisLabelFontFile = ''
      FTLE_FILEDisplay.PolarAxes.PolarAxisLabelBold = 0
      FTLE_FILEDisplay.PolarAxes.PolarAxisLabelItalic = 0
      FTLE_FILEDisplay.PolarAxes.PolarAxisLabelShadow = 0
      FTLE_FILEDisplay.PolarAxes.PolarAxisLabelFontSize = 12
      FTLE_FILEDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
      FTLE_FILEDisplay.PolarAxes.LastRadialAxisTextOpacity = 1.0
      FTLE_FILEDisplay.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
      FTLE_FILEDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
      FTLE_FILEDisplay.PolarAxes.LastRadialAxisTextBold = 0
      FTLE_FILEDisplay.PolarAxes.LastRadialAxisTextItalic = 0
      FTLE_FILEDisplay.PolarAxes.LastRadialAxisTextShadow = 0
      FTLE_FILEDisplay.PolarAxes.LastRadialAxisTextFontSize = 12
      FTLE_FILEDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
      FTLE_FILEDisplay.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
      FTLE_FILEDisplay.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
      FTLE_FILEDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''
      FTLE_FILEDisplay.PolarAxes.SecondaryRadialAxesTextBold = 0
      FTLE_FILEDisplay.PolarAxes.SecondaryRadialAxesTextItalic = 0
      FTLE_FILEDisplay.PolarAxes.SecondaryRadialAxesTextShadow = 0
      FTLE_FILEDisplay.PolarAxes.SecondaryRadialAxesTextFontSize = 12
      FTLE_FILEDisplay.PolarAxes.EnableDistanceLOD = 1
      FTLE_FILEDisplay.PolarAxes.DistanceLODThreshold = 0.7
      FTLE_FILEDisplay.PolarAxes.EnableViewAngleLOD = 1
      FTLE_FILEDisplay.PolarAxes.ViewAngleLODThreshold = 0.7
      FTLE_FILEDisplay.PolarAxes.SmallestVisiblePolarAngle = 0.5
      FTLE_FILEDisplay.PolarAxes.PolarTicksVisibility = 1
      FTLE_FILEDisplay.PolarAxes.ArcTicksOriginToPolarAxis = 1
      FTLE_FILEDisplay.PolarAxes.TickLocation = 'Both'
      FTLE_FILEDisplay.PolarAxes.AxisTickVisibility = 1
      FTLE_FILEDisplay.PolarAxes.AxisMinorTickVisibility = 0
      FTLE_FILEDisplay.PolarAxes.ArcTickVisibility = 1
      FTLE_FILEDisplay.PolarAxes.ArcMinorTickVisibility = 0
      FTLE_FILEDisplay.PolarAxes.DeltaAngleMajor = 10.0
      FTLE_FILEDisplay.PolarAxes.DeltaAngleMinor = 5.0
      FTLE_FILEDisplay.PolarAxes.PolarAxisMajorTickSize = 0.0
      FTLE_FILEDisplay.PolarAxes.PolarAxisTickRatioSize = 0.3
      FTLE_FILEDisplay.PolarAxes.PolarAxisMajorTickThickness = 1.0
      FTLE_FILEDisplay.PolarAxes.PolarAxisTickRatioThickness = 0.5
      FTLE_FILEDisplay.PolarAxes.LastRadialAxisMajorTickSize = 0.0
      FTLE_FILEDisplay.PolarAxes.LastRadialAxisTickRatioSize = 0.3
      FTLE_FILEDisplay.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
      FTLE_FILEDisplay.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
      FTLE_FILEDisplay.PolarAxes.ArcMajorTickSize = 0.0
      FTLE_FILEDisplay.PolarAxes.ArcTickRatioSize = 0.3
      FTLE_FILEDisplay.PolarAxes.ArcMajorTickThickness = 1.0
      FTLE_FILEDisplay.PolarAxes.ArcTickRatioThickness = 0.5
      FTLE_FILEDisplay.PolarAxes.Use2DMode = 0
      FTLE_FILEDisplay.PolarAxes.UseLogAxis = 0

      # reset view to fit data
      renderView1.ResetCamera()

      #changing interaction mode based on data extents
      renderView1.CameraPosition = [22.5, 3.0, 10000.0]
      renderView1.CameraFocalPoint = [22.5, 3.0, 0.0]

      # update the view to ensure updated data information
      renderView1.Update()

      # set scalar coloring
      ColorBy(FTLE_FILEDisplay, ('POINTS', 'FTLE - No T'))

      # rescale color and/or opacity maps used to include current data range
      FTLE_FILEDisplay.RescaleTransferFunctionToDataRange(True, False)

      # show color bar/color legend
      FTLE_FILEDisplay.SetScalarBarVisibility(renderView1, True)

      # get color transfer function/color map for 'FTLENoT'
      fTLENoTLUT = GetColorTransferFunction('FTLENoT')
      fTLENoTLUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
      fTLENoTLUT.InterpretValuesAsCategories = 0
      fTLENoTLUT.AnnotationsInitialized = 0
      fTLENoTLUT.ShowCategoricalColorsinDataRangeOnly = 0
      fTLENoTLUT.RescaleOnVisibilityChange = 0
      fTLENoTLUT.EnableOpacityMapping = 0
      fTLENoTLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 2.4434303303593876, 0.865003, 0.865003, 0.865003, 4.886860660718775, 0.705882, 0.0156863, 0.14902]
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
      fTLENoTPWF.Points = [0.0, 0.0, 0.5, 0.0, 4.886860660718775, 1.0, 0.5, 0.0]
      fTLENoTPWF.AllowDuplicateScalars = 1
      fTLENoTPWF.UseLogScale = 0
      fTLENoTPWF.ScalarRangeInitialized = 1

      # hide color bar/color legend
      FTLE_FILEDisplay.SetScalarBarVisibility(renderView1, False)

      # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
      fTLENoTLUT.ApplyPreset('jet', True)

      # Rescale transfer function
      fTLENoTLUT.RescaleTransferFunction(0.0, 5.0)

      # Rescale transfer function
      fTLENoTPWF.RescaleTransferFunction(0.0, 5.0)

      # create a new 'Legacy VTK Reader'
      CLOT_FILE = LegacyVTKReader(FileNames=[structFileName])

      # show data in view
      CLOT_FILEDisplay = Show(CLOT_FILE, renderView1)

      # trace defaults for the display properties.
      CLOT_FILEDisplay = Show(CLOT_FILE, renderView1)

      # trace defaults for the display properties.
      CLOT_FILEDisplay.Representation = 'Surface'
      CLOT_FILEDisplay.AmbientColor = [0.0, 0.0, 0.0]
      CLOT_FILEDisplay.ColorArrayName = [None, '']
      CLOT_FILEDisplay.OSPRayScaleArray = 'vertexNormals'
      CLOT_FILEDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      CLOT_FILEDisplay.SelectOrientationVectors = 'None'
      CLOT_FILEDisplay.ScaleFactor = 1.3261099815368653
      CLOT_FILEDisplay.SelectScaleArray = 'None'
      CLOT_FILEDisplay.GlyphType = 'Arrow'
      CLOT_FILEDisplay.GlyphTableIndexArray = 'None'
      CLOT_FILEDisplay.GaussianRadius = 0.06630549907684326
      CLOT_FILEDisplay.SetScaleArray = ['POINTS', 'vertexNormals']
      CLOT_FILEDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      CLOT_FILEDisplay.OpacityArray = ['POINTS', 'vertexNormals']
      CLOT_FILEDisplay.OpacityTransferFunction = 'PiecewiseFunction'
      CLOT_FILEDisplay.DataAxesGrid = 'GridAxesRepresentation'
      CLOT_FILEDisplay.PolarAxes = 'PolarAxesRepresentation'

      # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
      CLOT_FILEDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
      CLOT_FILEDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

      # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
      CLOT_FILEDisplay.DataAxesGrid.XTitleColor = [0.0, 0.0, 0.0]
      CLOT_FILEDisplay.DataAxesGrid.YTitleColor = [0.0, 0.0, 0.0]
      CLOT_FILEDisplay.DataAxesGrid.ZTitleColor = [0.0, 0.0, 0.0]
      CLOT_FILEDisplay.DataAxesGrid.GridColor = [0.0, 0.0, 0.0]
      CLOT_FILEDisplay.DataAxesGrid.XLabelColor = [0.0, 0.0, 0.0]
      CLOT_FILEDisplay.DataAxesGrid.YLabelColor = [0.0, 0.0, 0.0]
      CLOT_FILEDisplay.DataAxesGrid.ZLabelColor = [0.0, 0.0, 0.0]

      # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
      CLOT_FILEDisplay.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
      CLOT_FILEDisplay.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
      CLOT_FILEDisplay.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
      CLOT_FILEDisplay.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]

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

      #### saving camera placements for all active views

      # current camera placement for renderView1
      renderView1.InteractionMode = '2D'
      renderView1.CameraPosition = [22.5, 3.0, 10000.0]
      renderView1.CameraFocalPoint = [22.5, 3.0, 0.0]
      renderView1.CameraParallelScale = 12.813060868709753

      #### uncomment the following to render all views
      # RenderAllViews()
      # alternatively, if you want to write images, you can use SaveScreenshot(...).

      # destroy CLOT_FILE
      Delete(CLOT_FILE)
      del CLOT_FILE

      # set active source
      SetActiveSource(FTLE_FILE)

      # destroy FTLE_FILE
      Delete(FTLE_FILE)
      del FTLE_FILE

      #### saving camera placements for all active views

      # current camera placement for renderView1
      renderView1.InteractionMode = '2D'
      renderView1.CameraPosition = [22.5, 3.0, 10000.0]
      renderView1.CameraFocalPoint = [22.5, 3.0, 0.0]
      renderView1.CameraParallelScale = 12.813060868709753

      print('Wrote screenshot to %s' %outName)


