function poly_idx = polygon_approx(contour,dth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polygon approximation by fitting the curve with a sequence of lines
% contour: coordinates of contour pixels as an Nx2 matrix
% dth: threshold distance between the curve and the fitted line,
% smaller dth means finer approximation (default: 1)
% poly_idx: index of polygon pixels ( to use: contour(poly_idx,:) )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 1
    dth = 1;
end
contour = [contour; contour(1,:)];
N = size(contour,1);
poly_idx = zeros(N,1);
dth2 = dth^2;

poly_idx(1) = 1;
p1 = 1;
p_between = nan(N,1);
p_between(1) = p1 + 1;
len_between = 1;
p3 = p1 + 2;

while(p3<=N)
    p2 = p_between(1:len_between);
    % compute distance between p2 and line p1_p3
    if contour(p1,:) == contour(p3,:)
        d = (contour(p1,1)-contour(p2,1)).^2 + (contour(p1,2)-contour(p2,2)).^2;
        [dmax,imax] = max(d);
    else
        d = ((contour(p3,1)-contour(p2,1)).*(contour(p1,2)-contour(p2,2)) - ...
            (contour(p3,2)-contour(p2,2)).*(contour(p1,1)-contour(p2,1))).^2/...
            ((contour(p1,1)-contour(p3,1))^2+(contour(p1,2)-contour(p3,2))^2);
        [dmax,imax] = max(d);
    end
    % check if exceed threshold
    if dmax < dth2
        len_between = len_between + 1;
        p_between(len_between) = p3;
        p3 = p3 + 1;
    else
        p1 = p2(imax);
        p3 = p1+2;
        p_between = nan(N,1);
        p_between(1) = p1+1;
        len_between = 1;
        poly_idx(p1) = 1;
    end
end

poly_idx = logical(poly_idx);

end