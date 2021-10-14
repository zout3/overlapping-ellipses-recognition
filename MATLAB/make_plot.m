function out = make_plot(img, candidate)
img = dip_image(img,'bin');
%out = dip_array(img & ~berosion(img, 1, 1, 0));
out = dip_array(img);
x = dip_array(xx(img));
x = x - x(1,1) + 1;
y = dip_array(yy(img));
y = y - y(1,1) + 1;
for i = 1:size(candidate,1)
    cx = candidate(i, 1);
    cy = candidate(i, 2);
    a = candidate(i, 3);
    b = candidate(i, 4);
    d = candidate(i, 5);
    ellipse = dip_image((((x-cx)*cosd(d)+(y-cy)*sind(d))/a).^2 + ...
        (((y-cy)*cosd(d)-(x-cx)*sind(d))/b).^2 < 1);
    out(logical(dip_array(ellipse&~berosion(ellipse, 1, 1, 0)))) = 2;
end
out = uint8(dip_array(stretch(dip_image(out))));
end

