function [t_crossing] = detect_crossing(t,y, y_th, sign_th)

% Min-Young Kim, June 8, 2006
% To detect crossing points from below to above threshold in real vector y
% event triggering: y = y_th, with sign( y'(y_th) ) = sign_th 
% sign_th = 1 for positive slope, -1 for negative slope
% t_crossing is estimated using linear interpolation
% If (nx) is one of the crossing point, 
% t_crossing(nx) = t(nx-1) + (y_th-y(nx-1)/(y(nx) - y(nx-1))*delta_t

[n_row, n_col] = size(y);
if n_col == 1
    y = y'; t = t';
end;            % making a row vector

dy = sign_th*(y-y_th);
ind = find( ([0 dy]<0) & ([dy 0]>=0)); % indice right after crossing 
ind = ind(2:end-1);     % remove the first and the last crossing

delta_t= t(2)-t(1);
t_crossing = t(ind-1) + (y_th - y(ind-1))./(y(ind)- y(ind-1))*delta_t;
