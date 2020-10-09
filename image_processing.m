filenamestr = "../182804_e2_1um_10spf_2x63x_vary561_1pcto1p5pc.tif";
num_time = size(imfinfo(filenamestr),1);  % number of frames 
load('pair_positions.mat')
run_size = 2%num_time;

centrosome_pair_positions = cell(num_time, 1);


for i = 1:run_size
    clf;
    save_filename = strcat('labeled_points_peakfind/img_', sprintf('%03d',i), '.png');
    img = imread(filenamestr, i);
    orgimg = img;
    img = imgaussfilt(img, 3);    % Gaussian filter
    img = wiener2(img, [3 3]);    % Denoise
    p=FastPeakFind(img);
    %imagesc(img); hold on
    grayscaleIM(:,:,3) = orgimg;
    grayscaleIM(:,:,2) = orgimg;
    grayscaleIM(:,:,1) = orgimg;
    imagesc(grayscaleIM);
    hold on;
    axis equal
    axis off
    plot(p(1:2:end),p(2:2:end),'r+')
    
    if ~isempty(centrosome_pair_positions{i})
        list = centrosome_pair_positions{i};
    else
        list = [];
        list( :,1 ) = p(1:2:end);
        list(:, 2) = p(2:2:end);
    end
    
    
    roi = drawrectangle;
    tf = inROI(roi, list(:, 1), list(:, 2));
    
    
    
    %saveas(gcf,save_filename)
    %surf(1:size(img, 1), 1:size(img,2), img)
    
    
    centrosome_pair_positions{i} = list;
end


%% manual tracking
% % run_again = false;
% % 
% % i = 1;   
% % while i <= run_size
% %     clf;
% %     if ~isempty(centrosome_pair_positions{i})
% %         list = centrosome_pair_positions{i};   
% %     else
% %         list = [];
% %     end
% %     
% %   
% %    img = imread(filenamestr, i); 
% %    img = imgaussfilt(img, 2);    % Gaussian filter
% %    img = wiener2(img, [5 5]);    % Denoise
% %    
% %    grayscaleIM(:,:,3) = img;
% %    grayscaleIM(:,:,2) = img;
% %    grayscaleIM(:,:,1) = img;
% %    imagesc(grayscaleIM);
% %    axis equal
% %    axis off
% %    hold on
% %    if ~isempty(list)
% %        scatter(list(:,1), list(:,2),'b','filled'); % revise
% %        scatter(list(:,3), list(:,4),'r','filled');
% %    end
% %    
% %    if run_again 
% %       % generate list from k to run_size; 
% %       % for every point, generate a box where the point is at the center,
% %       % phase correlation
% %       for j = k: run_size-1
% %           imgj = imread(filenamestr, j);
% %           imgj = imgaussfilt(imgj, 2);    % Gaussian filter
% %           imgj = wiener2(imgj, [5 5]);    % Denoise
% %           imgjplus1 = imread(filenamestr, j+1);
% %           imgjplus1 = imgaussfilt(imgjplus1, 2);    % Gaussian filter
% %           imgjplus1 = wiener2(imgjplus1, [5 5]);    % Denoise
% %           new_list = zeros(size(list'));
% %           for jj = 1: 2: size(list,1)*4
% %               list = list';
% %               x = list(ind2sub(size(list), jj));
% %               y = list(ind2sub(size(list), jj+1));
% %               
% %               boxj_jj = imgj(round(x)-10:round(x)+10, round(y)-10:round(y)+10);
% %               boxjplus1_jj = imgjplus1(round(x)-10:round(x)+10, round(y)-10: round(y)+10);
% %               
% %               [shiftx, shifty] = xcorr2fft(boxj_jj, boxjplus1_jj);
% %               
% %               new_list(ind2sub(size(list),jj)) = x+shifty;
% %               new_list(ind2sub(size(list),jj+1)) = y+shiftx;
% %           end
% %           centrosome_pair_positions{j+1} = new_list';
% %       end
% %    end
% % %    
% % 
% %  
% %    waitforbuttonpress;
% %    value = double(get(gcf,'CurrentCharacter'));
% %    % press leftarrow to go backward, rightarrow to go forward. press up to
% %    % do something, right click 
% %    switch value
% %        case 28 % 28 leftarrow
% %            if i ~= 1
% %                i = i-1;
% %            end
% %        case 29 
% %            i = i+1;
% %        otherwise
% %            k = i;
% %            run_again = true;
% %            list = click_functions(value, list);                  
% %            centrosome_pair_positions{i} = list;
% %    end
% % end
% % 
% % %save("pair_positions.mat", 'centrosome_pair_positions');
% % 
% % function update_list = click_functions(value, list)
% %    % press leftarrow to go backward, rightarrow to go forward. press up to
% %    % do something, right click 
% %    switch value
% %        case 30 % 30 uparrow
% %            button = 3;
% %            while button ~= 1
% %                [x, y, button] = ginput(2);
% %                if button(1) == 3 && button(2) == 3 %right click twice
% %                    % add the two points to the cell position list
% %                    list = add_pair(x, y, list);
% %                    update_list = list;
% %                    if ~isempty(list)
% %                        scatter(list(:,1), list(:,2),'b','filled'); % revise
% %                        scatter(list(:,3), list(:,4),'r','filled');
% %                        hold on;
% %                    end
% %                else
% %                    update_list = list;
% %                end
% %            end
% %            
% %        case 31 % 31 downarrow
% %            button = 3;
% %            while button ~= 1
% %                [x, y, button] = ginput(1);
% %                if button(1) == 3  %right click once
% %                    % delete one cell from the position list
% %                    list = delete_single(x, y, list);
% %                    update_list = list;
% %                end
% %            end
% %        otherwise
% %    end
% % end
% % 
% % 
% % function update_list = delete_single(x,y, list) % x, y are numbers
% %     
% % 
% % end
% % 
% % function update_list = add_pair(x, y, list) % x, y are 2x1 list
% %     [n, ~] = size(list);
% %     update_list = zeros(n+1, 4);
% %     update_list(1:n,:) = list;
% %     update_list(n+1, 1) = x(1);
% %     update_list(n+1, 2) = y(1);
% %     update_list(n+1, 3) = x(2);
% %     update_list(n+1, 4) = y(2);
% % end

