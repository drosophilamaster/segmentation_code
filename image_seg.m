%%

filenamestr = "../20191017181230_JGFPHRFP_1um_10spf_2x63x_oiloil_4pc1pc_max_z_RFP.tif";
segmentationstr = "../20191017181230_JGFPHRFP_1um_10spf_2x63x_oiloil_4pc1pc_max_z_RFP_Simple Segmentation__.h5";

segmentation_data = h5read(segmentationstr, '/exported_data');

%% 
num_timepoints = size(imfinfo(filenamestr),1);

objects_positions = cell(num_timepoints,1);
i = 1;
while i <= num_timepoints-1
    img1 = imread(filenamestr, i);
    
    %gauss_img = imgaussfilt(img1,3);
    %bw1 = edge(gauss_img, 'Canny');
    seg = ( segmentation_data(1,:, :, i) == 1);
    seg = reshape(seg, [size(seg, 2), size(seg, 3)]);
    seg = seg';
    centroids = get_centroids(seg);
    c_positions = cell2mat(centroids);

    resize_scale = 5;
    img1 = imresize(img1, resize_scale);
    [xsize_rescale, ysize_rescale] = size(img1);
    img2 = draw_points(c_positions, size(img1));
    imshow(img1);
%     hold on;
%     h = imshow(img2);
%     set(h, 'AlphaData', 0.5);
%     button = 1;
%     while button ~= 2
%         [x, y, button] = ginput(1);
%         if button == 1 %left
%             c_positions = add_point(c_positions, size(img1) , [x,y]);
%             img2 = draw_points(c_positions, size(img1));
%             hold off;
%             imshow(img1);
%             hold on;
%             h = imshow(img2);
%             set(h, 'AlphaData', 0.5);
%         elseif button == 3 %right
%             c_positions = delete_point(c_positions, size(img1) , [x,y]);
%             img2 = draw_points(c_positions, size(img1));
%             hold off;
%             imshow(img1);
%             hold on;
%             h = imshow(img2);
%             set(h, 'AlphaData', 0.5);
%         end
%     end 
    c_positions = rmmissing(c_positions);
    objects_positions{i} = c_positions;
    hold off;
    i = i+1;
end


save('position_list_1.mat','objects_positions')


%%% functions
function updated_list = add_point(pos_list, img_size, click_position)
    [n, ~] = size(pos_list);
    updated_list = zeros(n+1, 2);
    sizeX = img_size(1);
	sizeY = img_size(2);
    new_x = click_position(1)/sizeX;
    new_y = click_position(2)/sizeY;
    updated_list(1:n,:) = pos_list;
    updated_list(n+1, :) = [new_x, new_y];
end

function updated_list = delete_point(pos_list, img_size, click_position)
    sizeX = img_size(1);
	sizeY = img_size(2);
    temp_positions(:,1) = pos_list(:,1)*sizeX;
    temp_positions(:,2) = pos_list(:,2)*sizeY;

    
    for k = 1:size(pos_list, 1)
        if abs(click_position(1)-temp_positions(k, 1)) <= 5 && abs(click_position(2) - temp_positions(k, 2)) <= 5
            pos_list(k,1) = NaN;
            pos_list(k,2) = NaN;
        end
    end
    
    updated_list = pos_list;
end


function img = draw_points(pos_list, img_size)
    sizeX = img_size(1);
	sizeY = img_size(2);
    img = zeros(sizeX, sizeY, 3);
    temp_positions(:,1) = pos_list(:,1)*sizeX;
    temp_positions(:,2) = pos_list(:,2)*sizeY;
    for k = 1:size(pos_list, 1)
        if isnan(temp_positions(k, 1))
        else
            position = round(temp_positions(k,:));
            pos_x= position(1);
            pos_y = position(2);
            img(pos_y-5: pos_y+5, pos_x-5:pos_x+5, 1) = 10;
        end
    end
end



%% get centriods from bw
function centroid_list = get_centroids(bw)
    bw_result = bw;
    x_size = size(bw, 1);
    y_size = size(bw, 2);
    centroid_list = {};

    bw = imfill(bw, 'holes');
    
    [B,L] = bwboundaries(bw, 'noholes');

    stats = regionprops(L,'Area','Centroid');
    a = zeros(length(B),1);
    for k = 1:length(B)
        boundary = B{k};
        area = stats(k).Area;
        a(k) = area-size(boundary, 1);
        % plot centroid if there is one
        if a(k) > 0
            centroid = stats(k).Centroid;
            centroid(1) = centroid(1)/x_size;
            centroid(2) = centroid(2)/y_size;
            %plot(centroid(1),centroid(2), 'ro');
            centroid_list{end+1,1} = centroid; 
        else
            bw_result(L == k) = 0;
        end
    end
end


