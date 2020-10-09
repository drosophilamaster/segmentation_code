filenamestr = "../182804_e2_1um_10spf_2x63x_vary561_1pcto1p5pc.tif";
num_time = size(imfinfo(filenamestr),1);  % number of frames 
load('pair_positions_2.mat');
run_size = 150 %num_time;

%centrosome_pair_positions = cell(num_time, 1);

i = 1;
run_again = false;
while i <= run_size
    clf;
    save_filename = strcat('peakfind_fixed/img_', sprintf('%03d',i), '.png');
    img = imread(filenamestr, i);
    orgimg = img;
    img = imgaussfilt(img, 3);    % Gaussian filter
    img = wiener2(img, [3 3]);    % Denoise
   % p=FastPeakFind(img);
    %imagesc(img); hold on
    grayscaleIM(:,:,3) = orgimg;
    grayscaleIM(:,:,2) = orgimg;
    grayscaleIM(:,:,1) = orgimg;
    imagesc(grayscaleIM);
    hold on;
    axis equal
    axis off
    text(0,0, num2str(i), 'Color', 'r')
    
    if ~isempty(centrosome_pair_positions{i})
        list = centrosome_pair_positions{i};
    else
        list = [];
%         list( :,1 ) = p(1:2:end);
%         list(:, 2) = p(2:2:end);
    end
    scatter(list(:,1),list(:,2),'r', 'filled');
    
%     roi = drawrectangle;
%     tf = inROI(roi, list(:, 1), list(:, 2));
%     
    


waitforbuttonpress;
   value = double(get(gcf,'CurrentCharacter'));
   % press leftarrow to go backward, rightarrow to go forward. press up to
   % do something, right click 
   switch value
       case 28 % 28 leftarrow
           if i ~= 1
               i = i-1;
           end
       case 29 
           i = i+1;
       case 30
           k = i;
           run_again = true;
           roi = drawrectangle;
           tf = inROI(roi, list(:, 1), list(:, 2));
           list = delete_points(list, tf);
           centrosome_pair_positions{i} = list;
           
       otherwise
           k = i;
           run_again = true;
           [x, y, button] = ginput(1)
           list = add_points(list, x,y);
           centrosome_pair_positions{i} = list;
   end
    
    %saveas(gcf,save_filename)
    %surf(1:size(img, 1), 1:size(img,2), img)
    
    %i = i+1;
end

save("pair_positions_2.mat", 'centrosome_pair_positions');



function new_list = delete_points(list, tf)
    list(tf, :) = nan();
    new_list = rmmissing(list);
end


function new_list = add_points(list, x, y)
    [n, ~] = size(list);
    new_list = zeros(n+1, 2);
    new_list(1:n,:) = list;
    new_list(n+1, 1) = x
    new_list(n+1, 2) = y
end