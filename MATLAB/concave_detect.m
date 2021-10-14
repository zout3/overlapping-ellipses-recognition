function concave_idx = concave_detect(poly, type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concave points detection on ploygon approximation by checking the angle
% btween every two connected sides
% poly: coordinates of polygon pixels as an Nx2 matrix
% type: 1 is an outer polygon, -1 is an inner polygon
% concave_idx: index of concave pixels ( to use: poly(concave_idx,:) )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1
    type = 1;
end
N = size(poly,1);
concave_idx = zeros(1,N);

v1 = [poly(1,1)-poly(N,1);poly(1,2)-poly(N,2)];
v2 = [poly(2,1)-poly(1,1);poly(2,2)-poly(1,2)];
if type*(v1(1)*v2(2) - v1(2)*v2(1)) < 0
    concave_idx(1) = 1;
end

for i = 1:(N-2)
    v1 = [poly(i+1,1)-poly(i,1);poly(i+1,2)-poly(i,2)];
    v2 = [poly(i+2,1)-poly(i+1,1);poly(i+2,2)-poly(i+1,2)];
    if type*(v1(1)*v2(2) - v1(2)*v2(1)) < 0
        concave_idx(i+1) = 1;
    end
end

v1 = [poly(N,1)-poly(N-1,1);poly(N,2)-poly(N-1,2)];
v2 = [poly(1,1)-poly(N,1);poly(1,2)-poly(N,2)];
if type*(v1(1)*v2(2) - v1(2)*v2(1)) < 0
    concave_idx(N) = 1;
end

concave_idx = logical(concave_idx);

end




