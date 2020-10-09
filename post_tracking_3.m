addpath('../')
load('new_index_track.mat'); % tracks, 4 column [x,y,time,index], sorted in index
filenamestr = "../182804_e2_1um_10spf_2x63x_vary561_1pcto1p5pc.tif";
load("positionshift.mat");
start_time = 1;
run_time = 63;


xshift = cumsum([0;positionshift(:,1)]);
yshift = cumsum([0;positionshift(:,2)]);

maximum_index = max(tracks_new(:,5));
centrosome_pairs = [1:2:maximum_index; 2:2:maximum_index]';
distance = nan(run_time-start_time, size(centrosome_pairs,1));

for i = start_time: run_time
     ilist = tracks_new((tracks_new(:,3) == i),:);
    for k = 1:size(centrosome_pairs, 1)
        
        x1 = ilist(ilist(:,5) == centrosome_pairs(k, 1), 1);
        x2 = ilist(ilist(:,5) == centrosome_pairs(k, 2), 1);
        y1 = ilist(ilist(:,5) == centrosome_pairs(k, 1), 2);
        y2 = ilist(ilist(:,5) == centrosome_pairs(k, 2), 2);
        if ~isempty(x1) && ~isempty(x2) && ~isempty(y1) && ~isempty(y1)
            distance(i, k) = sqrt((x1-x2)^2+(y1-y2)^2)*0.361984;
        end
    end
end
    %%
plot(start_time:10:run_time*10, nanmean(distance,2), 'color', 'r', 'LineWidth', 5)
hold on;
plot(start_time:10:run_time*10, distance, '+')

xlabel('time [s]')
ylabel('distance [micron]')