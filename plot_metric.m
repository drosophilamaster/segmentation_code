
cutoff =0.87;

for i = -pi/2:0.01:0
    clf
    viscircles([0,0],1,'Color','g');
    hold on
    viscircles([0,0],cutoff,'Color','g');
    y = sin(i);
    x1 = -cos(i);
    x2 = cos(i);

    plot( x1, y,'or','MarkerSize',5,'MarkerFaceColor','r');  % point1
    plot(x2,y,'or','MarkerSize',5,'MarkerFaceColor','r');   % point2
    plot([x1 x2], [y y])    % side 12
    
    
    % insert point3
    point3 = compute_point3(y, abs(x1), cutoff);
    x3 = point3(1);
    y3 = point3(2);
    
    
    plot(x3,y3,'or','MarkerSize',5,'MarkerFaceColor','r');
    % plot side 23, 13
    plot([x1 x3], [y y3])
    plot([x2 x3], [y y3])
    
    x_m = mean([x1, x2, x3]);
    y_m = mean([y, y, y3]);
    
    plot(x_m,y_m,'or','MarkerSize',5,'MarkerFaceColor','b');
    
    axis equal
   
    xlim([-1.5, 1.5])
    ylim([-1.5, 1.5])
    title(sprintf("|centriod-circumcenter|/circumradius: %f", cutoff))
    
    hold off;
    pause(0.001)
end


%%

function thirdpoint = compute_point3(y_value, x_value, cutoff)
    center_y = abs(y_value)*1/3+y_value;
    radius_small = sqrt((x_value/3)^2 + (center_y-y_value)^2); 
    intersection = gettwocircle(y_value, center_y, radius_small, cutoff);
    if isnan(intersection)
        thirdpoint = nan(1,2);
    else
        thirdpoint = get3rdpoint(y_value, intersection);
    end
end



function two_circle = gettwocircle(y_pos,center_y, r, cutoff) % return the intersection of two circles, take the left one
    %viscircles([0, center_y], r);
    syms x y;
    eqns = [x^2+y^2 == cutoff^2, x^2+(y-center_y)^2 == r^2];
    vars = [x,y];
    [cx, cy] = vpasolve(eqns, vars);
    if imag(cx(1)) ~= 0 || imag(cy(1))~= 0
        two_circle = nan(1,2);
    elseif cy(1) - y_pos >0 && cx(1) < 0
        two_circle = [cx(1), cy(1)];
    elseif cy(2) - y_pos >0 && cx(2) < 0
        two_circle = [cx(2), cy(2)];
    elseif cy(1) - y_pos <0 && cx(1) > 0
        two_circle = [cx(1), cy(1)];
    elseif cy(2) - y_pos <0 && cx(2) > 0
        two_circle = [cx(2), cy(2)];
    end
end




function third_point = get3rdpoint(y_pos, point) % return the position of the third point
    x0 = point(1);
    y0 = point(2);
    syms x y;
    eqns = [x^2+y^2 == 1, y-y_pos == x*(y0-y_pos)/x0];
    vars = [x,y];
    [sx, sy] = vpasolve(eqns, vars);
    if sy(1) - y_pos >0 && sx(1) < 0 && y0>y_pos
        third_point = [sx(1), sy(1)];
    elseif sy(2) - y_pos >0 && sx(2) < 0 && y0>y_pos
        third_point = [sx(2), sy(2)];
    elseif sy(1) - y_pos <0 && sx(1) > 0 && y0<y_pos
        third_point = [sx(1), sy(1)];
    elseif sy(2) - y_pos <0 && sx(2) > 0 && y0<y_pos
        third_point = [sx(2), sy(2)];
    else
        third_point = nan(1,2);
    end
end