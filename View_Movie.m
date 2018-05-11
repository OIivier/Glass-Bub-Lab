function rotorboolean = View_Movie(DataName);
[SF,Lx,Ly] = daRead12(DataName, 0);
% Data = SF(:,:,1:round(size(SF,3)/2));
Data = SF;
Data = bandpass(Data);
movie_maker_mod(Data,4000,4800,1,0,0,0,0,1,0,1);