%
%
% Thomas Quail
% Detrend the time series from Alex data
%
% 2015
load('nevinTS.txt');
%
%
ts=nevinTS;
slidingwindow=5000;
numbIts=floor(length(nevinTS)/slidingwindow);
%
StartPoints=linspace(slidingwindow,...
    slidingwindow*floor(length(nevinTS)/slidingwindow),floor(length(nevinTS)/slidingwindow));

% 
DetrendedTimeSeries=[];
for k=1:numbIts-1;
    CurrentVector=detrend(ts((k*slidingwindow+1):(k*slidingwindow+slidingwindow)));
    DetrendedTimeSeries=[DetrendedTimeSeries CurrentVector'];
end
