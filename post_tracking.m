addpath('../')
load('tracked_particles.mat'); % tracks, 4 column [x,y,time,index], sorted in index
filenamestr = "../182804_e2_1um_10spf_2x63x_vary561_1pcto1p5pc.tif";
load("positionshift.mat");
start_time = 1;
run_time = 63;

xshift = cumsum([0;positionshift(:,1)]);
yshift = cumsum([0;positionshift(:,2)]);
tracks_new = [tracks, zeros(size(tracks, 1), 1)];

i = start_time ;
while i <= run_time
    clf;
    
    clist = tracks_new((tracks_new(:,3) == i),:);
    img1 = imread(filenamestr , i);
    
    grayscaleIM(:,:,3) = img1;
    grayscaleIM(:,:,2) = img1;
    grayscaleIM(:,:,1) = img1;
    imagesc(grayscaleIM)
    axis equal;
    hold on;
    
    scatter((clist(:, 1)+yshift(i)), (clist(:,2)+xshift(i)), 'b','filled')
    text(clist(:, 1)+yshift(i), clist(:,2)+xshift(i), num2str(clist(:, 5)), 'Color', 'r')
    waitforbuttonpress;
    value = double(get(gcf,'CurrentCharacter'));
    switch value
        case 29 % forward, done
            i=i+1;
        case 30 % upward, do it again
            roi = drawfreehand;
            tf = inROI(roi, clist(:,1)+yshift(i), clist(:,2)+xshift(i));
            indeces = clist(tf, 4);
            if length(indeces) == 2
                tracks_new(tracks_new(:, 4) == indeces(1) & tracks_new(:, 3)>=i, 5) = max(tracks_new(:, 5))+1;
                tracks_new(tracks_new(:, 4) == indeces(2) & tracks_new(:, 3)>=i, 5) = max(tracks_new(:, 5))+1;
            end
        case 28 % backward,
            if i >= 1
                i=i-1;
            end
        otherwise
    end
end

save('new_index_track.mat','tracks_new');
