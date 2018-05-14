function GBRColorMap = GBRColorMap(ColorMapResolution)

GBRColorMapTop = flipud([zeros(1,ColorMapResolution);2.^[(1.0/ColorMapResolution):(1.0/ColorMapResolution):1]-2^(1.0/ColorMapResolution);zeros(1,ColorMapResolution)]');
GBRColorMapBottom = [2.^[(1.0/ColorMapResolution):(1.0/ColorMapResolution):1]-2^(1.0/ColorMapResolution);zeros(1,ColorMapResolution);zeros(1,ColorMapResolution)]';
GBRColorMapTop = flipud([zeros(1,ColorMapResolution);[(1.0/ColorMapResolution):(1.0/ColorMapResolution):1];zeros(1,ColorMapResolution)]');
GBRColorMapBottom = [[(1.0/ColorMapResolution):(1.0/ColorMapResolution):1];zeros(1,ColorMapResolution);zeros(1,ColorMapResolution)]';

GBRColorMap = [GBRColorMapTop;GBRColorMapBottom];