function elist = overlap(img, candidate, prcntl, cover_rate)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overlap candidate ellipses on the image in the order of their percentile
% score, and return ellipses on the top of the overlapping
% img: image as a 2-d matrix
% candidate: list of candidate ellipses (center coordinates, length, width
% and angle in degree) as an Nx5 matrix
% prcntl: percentiles used to compute ellipse scores (default:[30,50,70])
% elist: list of ellipses on the top of the overlapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

img = dip_image(img,'bin');
switch nargin
    case 2
        prcntl = [30,50,70];
        cover_rate = 0.95;
    case 3
        cover_rate = 0.95;
end


% distance transform of the contour, used to compute score
EDT = dip_array(dt(~(img&~berosion(img,1,1,0))));
X = dip_array(xx(img));
X = X - X(1,1) + 1;
Y = dip_array(yy(img));
Y = Y - Y(1,1) + 1;
[len_y, len_x] = size(X);
num_ellipse = size(candidate,1);
num_p = numel(prcntl);
IDX3D = len_y*len_x*(0:(num_p-1));

% the following computes a window for each ellipse
cx = candidate(:,1);
cy = candidate(:,2);
a = candidate(:,3);
b = candidate(:,4);
d = candidate(:,5);
dx = abs((a.^2.*cosd(d)+b.^2.*tand(d).*sind(d))./sqrt(a.^2+b.^2.*tand(d).^2));
dy = abs((b.^2.*cosd(d)+a.^2.*tand(d).*sind(d))./sqrt(b.^2+a.^2.*tand(d).^2));
dx(isnan(dx)) = b(isnan(dx));
dy(isnan(dy)) = a(isnan(dy));
x_min = max(floor(cx-dx),1);
y_min = max(floor(cy-dy),1);
x_max = min(ceil(cx+dx),len_x);
y_max = min(ceil(cy+dy),len_y);

% store score and area of ellipses on the overlap for different percentiles
scoreAll = Inf(len_y, len_x, num_p);
areaAll = zeros(len_y, len_x, num_p);
overlap = -1*ones(len_y, len_x, num_p);

for i = 1:num_ellipse
    % use the window to plot ellipse
    edt = EDT(y_min(i):y_max(i),x_min(i):x_max(i));
    x = X(y_min(i):y_max(i),x_min(i):x_max(i));
    y = Y(y_min(i):y_max(i),x_min(i):x_max(i));
    cx = candidate(i,1);
    cy = candidate(i,2);
    a = candidate(i,3);
    b = candidate(i,4);
    d = candidate(i,5);
    
    % use in and out to approximate ellipse contour
    in = (((x-cx)*cosd(d)+(y-cy)*sind(d))/a).^2 + ...
        (((y-cy)*cosd(d)-(x-cx)*sind(d))/b).^2 < 1;
    out = (((x-cx)*cosd(d)+(y-cy)*sind(d))/(a-1)).^2 + ...
        (((y-cy)*cosd(d)-(x-cx)*sind(d))/(b-1)).^2 > 1;
    score = prctile(edt(in & out), prcntl);
    
    % put windowed ellipse on image of original size
    elps = zeros(len_y, len_x);
    elps(y_min(i):y_max(i),x_min(i):x_max(i)) = in;
    idx = find(elps) + IDX3D;
    score =  repmat(score, [size(idx,1),1]);
    
    % overlap if new ellipse has lower score or same score but larger area
    cover_idx = (scoreAll(idx) > score) | ...
        ( (scoreAll(idx) == score) & (areaAll(idx) < a*b) );
    idx = idx(cover_idx);
    scoreAll(idx) = score(cover_idx);
    areaAll(idx) = a*b;
    overlap(idx) = i;
end

remain = zeros(num_p, num_ellipse);
minArea = (1-cover_rate)*pi*candidate(:,3).*candidate(:,4);
for i = 1:num_p
    measure = regionprops(overlap(:,:,i), 'Area');
    % get rid of ellipse with small covering area
    remain(i,1:numel(measure)) = cat(1,measure.Area) > minArea(1:numel(measure));
end
ID = repmat(1:num_ellipse,[num_p,1]);
elist = candidate(unique(ID(logical(remain))),:);
end

