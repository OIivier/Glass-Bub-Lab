%% SpectrHomemade Function
function Heatmap2(Data,Scale,Ticks)
hold off;
hm = heatmap(Data);
caxis(hm, [-1,1]*Scale);

if Ticks==0
    old_warning_state = warning('off', 'MATLAB:structOnObject');
    hs = struct(hm);
    warning(old_warning_state);
    %now we can get at the private properties XAxis and YAxis
    hs.XAxis.TickValues = [];
    hs.YAxis.TickValues = [];
end


hm.Colormap = GBRColorMap(1024);
