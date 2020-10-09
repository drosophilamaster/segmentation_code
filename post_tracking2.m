addpath('../')
load('new_index_track.mat'); % tracks, 4 column [x,y,time,index], sorted in index
filenamestr = "../182804_e2_1um_10spf_2x63x_vary561_1pcto1p5pc.tif";
load("positionshift.mat");
start_time = 1;
run_time = 63;

xshift = cumsum([0;positionshift(:,1)]);
yshift = cumsum([0;positionshift(:,2)]);


%%
for i = start_time: run_time
    save_filename = strcat('track_plot_pair/img_', sprintf('%03d',i), '.png');
    clf;
    img1 = imread(filenamestr, i);
    indices = find(tracks_new(:,3) == i);
    
    
    grayscaleIM(:,:,3) = img1;
    grayscaleIM(:,:,2) = img1;
    grayscaleIM(:,:,1) = img1;
    imagesc(grayscaleIM)
    axis equal;
    axis off;
    hold on; 
    
    scatter((tracks_new(indices, 1)+yshift(i)), (tracks_new(indices,2)+xshift(i)), 'r','+')
    
    ilist = tracks_new((tracks_new(:,3) == i),:);
    maximum_index = max(ilist(:,5));
    centrosome_pairs = [1:2:maximum_index; 2:2:maximum_index]';
    
    for k = 1:size(centrosome_pairs, 1)
        x1 = ilist(ilist(:,5) == centrosome_pairs(k, 1), 1);
        x2 = ilist(ilist(:,5) == centrosome_pairs(k, 2), 1);
        y1 = ilist(ilist(:,5) == centrosome_pairs(k, 1), 2);
        y2 = ilist(ilist(:,5) == centrosome_pairs(k, 2), 2);
        
        plot([x1,x2]+yshift(i), [y1,y2]+xshift(i), 'g')
    end
    %saveas(gcf, save_filename);
    pause(0.2)
end