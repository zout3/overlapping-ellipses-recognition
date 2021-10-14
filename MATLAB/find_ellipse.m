function candidate = find_ellipse(img, ratio, angle, minArea, nMax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract possible ellipses using compression and Euclidian distance transform
% img: image as a 2-d matrix
% ratio: searching grid for ratio of length over width (default: 1:0.2:3)
% angle: searching grid for ellipse angle in degree (default: 5:5:180)
% minArea: ellipses with area under minArea will be dropped (default: sum(img)/50) 
% nMax: maximum number of candidate ellipses (default: 1e5)
% candidate: list of candidate ellipses (center coordinates, length, width
% and angle in degree) as an Nx5 matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
img = dip_image(img,'bin');
switch nargin
    case 1
        ratio = 1:0.2:4;
        angle = 5:5:180;
        nMax = 1e6;
        minArea = sum(img)/50;
    case 2
        angle = 5:5:180;
        nMax = 1e6;
        minArea = sum(img)/50;
    case 3
        nMax = 1e6;
        minArea = sum(img)/50;
    case 4
        nMax = 1e6;
end
candidate = zeros(nMax,5);
count = 0;

x_centered = dip_array(xx(img));
y_centered = dip_array(yy(img));
ref = [x_centered(1,1)-1; y_centered(1,1)-1];
[len_y, len_x] = size(x_centered);
coordinate = [x_centered(img)';y_centered(img)']; % pixel coordinates with
                                                  % origin at the center  
num_pixel = size(coordinate, 2);
num_angle = numel(angle);
num_ratio = numel(ratio);
num_pair = num_angle * num_ratio;
angle = reshape(repmat(angle',[1,num_ratio]),[1,num_pair]);
ratio = reshape(repmat(ratio,[num_angle,1]),[1,num_pair]);

% construct compressing transform matrics with different parameters
cos_angle = cosd(angle);
sin_angle = sind(angle);
cos2 = cos_angle.^2;
sin2 = sin_angle.^2;
cossin = cos_angle.*sin_angle;
Ta = (1./ratio - 1);
Tinva = (ratio - 1);
T = zeros(2,2,num_pair);
T(1,1,:) = Ta.*cos2 +1;
T(1,2,:) = Ta.*cossin;
T(2,1,:) = T(1,2,:);
T(2,2,:) = Ta.*sin2 +1;
Tinv = zeros(2,2,num_pair);
Tinv(1,1,:) = Tinva.*cos2 +1;
Tinv(1,2,:) = Tinva.*cossin;
Tinv(2,1,:) = Tinv(1,2,:);
Tinv(2,2,:) = Tinva.*sin2 +1;

% conduct compressing transform
T_coordinate = zeros(2,num_pixel,num_pair);
for i = 1:num_pair
    T_coordinate(:,:,i) = round(T(:,:,i)*coordinate);
end
T_coordinate = T_coordinate - ref; % with origin at top left

% modify coordinates due to image size change by compressing
up_left = squeeze(min(T_coordinate,[],2)) - [1; 1];
low_right = squeeze(max(T_coordinate,[],2)) - [len_x;len_y];
T_coordinate = T_coordinate - ...
    permute(repmat(up_left,[1,1,num_pixel]),[1,3,2]);
image_size = low_right - up_left + [len_x; len_y];

% drop parameters with 1-d image after compressing
id_param = 1:num_pair;
id_param = id_param(min(image_size)~=1);

for i = id_param
    % construct compressed contour image
    size_i = [image_size(2,i),image_size(1,i)];
    edt = zeros(size_i);
    edt(sub2ind(size_i,T_coordinate(2,:,i),T_coordinate(1,:,i))) = 1;
    % conduct EDT and find local maximum
    edt = dt(dip_image(edt,'bin'),0);
    % smooth the EDT to extract local maximum
    smth = unif(edt, 3);
    localmax = find(maxima(smth, 2));
    [y, x] = ind2sub(size_i,localmax + 1);
    n = numel(x);
    if count + n > nMax
        error('Error: candidate ellipse number exceeds nMax');
    end
    % acquire original coordinates by the inverse of compressing
    cord_new = Tinv(:,:,i)*([x';y'] + up_left(:,i) + ref) - ref;
    candidate((count+1):(count+n),[1,2]) = cord_new';
    candidate((count+1):(count+n),3) = ratio(i);
    candidate((count+1):(count+n),4) = double(edt(localmax));
    candidate((count+1):(count+n),5) = angle(i);
    count = count + n;
end
candidate(:,[1,2]) = round(candidate(:,[1,2]));
candidate(:,3) = candidate(:,3).*candidate(:,4);
candidate = candidate(1:count,:);

area = pi*candidate(:,3).*candidate(:,4);
candidate = candidate(area>minArea,:);
end






