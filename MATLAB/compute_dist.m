function distMat = compute_dist(img, elist, contour, seg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the distance matrix between the P contour segments and Q ellipses
% img: image as a 2-d matrix
% elist: list of ellipses (center coordinates, length, width
% and angle in degree) as an Nx5 matrix
% contour: structure containing coordinates of contour pixels
% seg: structure containing index of contour segment points
% distMat: PxQ distance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

distMat = [];
x = dip_array(xx(img));
x = x - x(1,1) + 1;
y = dip_array(yy(img));
y = y - y(1,1) + 1;
num_ellipse = size(elist,1);

for k = 1:numel(contour)
    num_pixel = size(contour(k).cord,1);
    dist = zeros(num_pixel, num_ellipse);
    for i = 1:num_ellipse
        cx = elist(i, 1);
        cy = elist(i, 2);
        a = elist(i, 3);
        b = elist(i, 4);
        d = elist(i, 5);
        elps = (((x-cx)*cosd(d)+(y-cy)*sind(d))/a).^2 +...
            (((y-cy)*cosd(d)-(x-cx)*sind(d))/b).^2 < 1;
        temp = dip_array(dt(~(elps & ~berosion(elps, 1, 1, 0))));
        dist(:,i) = temp(sub2ind(size(temp),contour(k).cord(:,2),contour(k).cord(:,1)));
    end
    
    num_segment = numel(seg(k).id);
    if num_segment < 2
        distMat = [distMat; sum(dist)];
        continue;
    end
    
    distMat_temp = zeros(num_segment, num_ellipse);
    for i = 1:(num_segment-1)
        distMat_temp(i,:) = sum(dist(seg(k).id(i):(seg(k).id(i+1)-1),:),1);
    end
    distMat_temp(num_segment,:) = sum(dist(seg(k).id(num_segment):end,:),1) + ...
        sum(dist(1:(seg(k).id(1)-1),:),1);
    distMat = [distMat; distMat_temp];
    
end

end

