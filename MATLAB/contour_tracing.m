function contour = contour_tracing(img, nMax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract contour sequence based on radial sweep method
% img: image as a 2-d matrix
% nMax: maximum number of contour pixels (default: 1e5)
% contour: structure containing coordinates of contour pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1
    nMax = 1e5;
end

% Add blank frame
img = [zeros(size(img,1),1),img,zeros(size(img,1),1)];
img = [zeros(1,size(img,2));img;zeros(1,size(img,2))];
img = dip_array(img&~berosion(logical(img), 1, 1, 0));

num_loop = 0;
while 1
    % Initial pixel
    [y1,x1] = find(img, 1, 'first');
    if isempty(y1)
        break
    end
    contour_list = nan(nMax,2);
    contour_list(1,:) = [x1,y1];
    count = 1;
    back = 1;
    
    % Neighbor index
    dx0 = [-1,0,1,1,1,0,-1,-1];
    dy0 = [-1,-1,-1,0,1,1,1,0];
    
    % Check if it's single pixel
    if isempty(find(img(sub2ind(size(img),y1+dy0,x1+dx0)),1))
        img(sub2ind(size(img),y1,x1)) = 0;
        continue;
    end
    
    while 1
        % Radial sweep
        dx = [dx0(back:end), dx0(1:(back-1))];
        dy = [dy0(back:end), dy0(1:(back-1))];
        x = contour_list(count,1) + dx;
        y = contour_list(count,2) + dy;
        idx = find(img(sub2ind(size(img),y,x)), 1, 'first');
        count = count + 1;
        if count > nMax
            error('Error: contour pixel number exceeds nMax');
        end
        contour_list(count,:) = [x(idx), y(idx)];
        back = mod(mod(back+idx-2,8)+5,8)+1;
        % Check stopping criterion
        if contour_list(1,:) == [x(idx), y(idx)]
            dx = [dx0(back:end), dx0(1:(back-1))];
            dy = [dy0(back:end), dy0(1:(back-1))];
            x = contour_list(count,1) + dx;
            y = contour_list(count,2) + dy;
            idx = find(img(sub2ind(size(img),y,x)), 1, 'first');
            if contour_list(2,:) == [x(idx),y(idx)]
                break;
            end
        end
    end
    contour_list = contour_list(1:(count-1),:);
    img(sub2ind(size(img),contour_list(:,2),contour_list(:,1))) = 0;
    if count > 10
        contour_list = contour_list - 1;
        num_loop = num_loop + 1;
        contour(num_loop).cord = contour_list;
    end
end

end




