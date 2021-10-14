clear;clc
run('D:\Program Files\DIPimage 2.9\dipstart.m');
image = imbinarize(rgb2gray(imread('pics/example.jpg')));
stats = regionprops(image, 'Image', 'BoundingBox');

dt = 1;  % distance threshold
lambda = 0.5;  % regularization parameter

elist = [];
for i = 1:length(stats)
    img = int8(stats(i).Image);

    % Extract contour coordinates (both outer and inner contours)
    contour = contour_tracing(img);

    seg = [];
    poly = [];
    for j = 1:numel(contour) % For each contour

        % Find polygon approximation points
        poly_idx = polygon_approx(contour(j).cord, dt);

        % Find concave points
        concave_idx = concave_detect(contour(j).cord(poly_idx,:), ~(j-1)*2-1);

        temp = 1:size(contour(j).cord,1);
        poly(j).id = temp(poly_idx);
        seg(j).id = poly(j).id(concave_idx);
    end

    % List all possible ellipses using distance transform
    candidate = find_ellipse(img);

    % Reduce the number of candidate ellipses using overlaying method
    candidate = overlap(img, candidate);

    % Compute distance between candidate ellipses and contour segments
    distMat = compute_dist(img, candidate, contour, seg);

    % Find optimal ellipses using integer programming
    sqrtArea = sqrt(sum(img(:)));
    ind = integer_progamming(distMat, lambda*sqrtArea);

    % Output
    imwrite(make_plot(img, candidate(ind, :)), ['pics/', int2str(i), '.jpg']);
    candidate(ind, 1) = candidate(ind, 1) + floor(stats(i).BoundingBox(1));
    candidate(ind, 2) = candidate(ind, 2) + floor(stats(i).BoundingBox(2));
    elist = [elist; candidate(ind, :)];
end

imwrite(make_plot(image, elist), 'pics/result.jpg');

