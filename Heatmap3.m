function Heatmap3(Data,ColorMap,ScaleMin,ScaleMax,ColorBar)

if ~exist('ColorMap','var')
	ColorMap = jet;
end

if ~exist('ScaleMin','var')
	ScaleMin = 0;
end

if ~exist('ScaleMax', 'var')
    ScaleMax = 1;
end

if ~exist('ColorBar', 'var')
    ColorBar = 1;
end

imagesc(Data,[ScaleMin,ScaleMax]);
colormap(gca,ColorMap);
if ColorBar==1
   colorbar;  
end