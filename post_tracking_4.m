
%%
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



% find the velocity of each particle.
x = tracks_new(:, 1);
y = tracks_new(:, 2);
difx = x(3:1:end)-x(1:1:end-2); difx = [nan(); difx; nan()];
dify = y(3:1:end)-y(1:1:end-2); dify = [nan(); dify; nan()];
difx(diff(tracks_new(:, 4)) == 1) = nan();
difx(diff([1;tracks_new(:, 4)]) == 1) = nan();
dify(diff(tracks_new(:, 4)) == 1) = nan();
dify(diff([1;tracks_new(:, 4)]) == 1) = nan();

tracks_new(:, 5) = difx;
tracks_new(:, 6) = dify;




a = optimvar('a');
b = optimvar('b');
xi_1 = optimvar('xi_1');
xi_2 = optimvar('xi_2');

%% loop to get the residue
for i = 1%start_time: run_time-1
    
%     
%     clf;
%     img1 = imread(filenamestr, i);
%     
%     grayscaleIM(:,:,3) = img1;
%     grayscaleIM(:,:,2) = img1;
%     grayscaleIM(:,:,1) = img1;
%     imagesc(grayscaleIM)
%     axis equal;
%     axis off
%     hold on; 
%  
    ilist = tracks_new((tracks_new(:,3) == i),:);
    P = ilist(:,1:2)+[yshift(i), xshift(i)];
    DT = delaunay(P);
    triplot(DT,P(:,1), P(:,2));
    
    
    TR = triangulation(DT,P);
    % maybe clean that a bit, to make it 
    
    
    conn_matrix = zeros(size(ilist, 1));
    pair_matrix = zeros(size(ilist, 1));
    %% looping the triagulation to get connectivity matrix
    for k = 1: size(ilist,1)
        for kk = 1: k-1
            if isConnected(TR, k, kk) % connectivity
                conn_matrix(k,kk) = 1;
            end
            if mod(k,2) == 0 && kk == k-1
                pair_matrix(k, kk) = 1;
            end
        end
    end

    conn_matrix = conn_matrix' + conn_matrix; % if conn_matrix (i,j) = 1, i and j are neighbors.
    pair_matrix = pair_matrix' + pair_matrix; % same but not neibors but pairs. 
   neighbor_matrix =  conn_matrix-pair_matrix ;
  
    %% compute the residue, a function of (a, b, c, d, x, y, vx, vy)
    
    for ii = 1: size(ilist,1) % looping all the particles
        if 1 ==1 % not boundary and velocity exist.
            for iii = 1: ii-1
                if neighbor_matrix(ii, iii) == 1
                    r = ilist(iii, 1:2) - ilist(ii, 1:2);
                elseif pair_matrix(ii, iii) == 1
                    r_p = ilist(iii, 1:2) - ilist(ii, 1:2);
                end
            end
        end
        % somegiganticfunction x,y  = somegiganticfunction +
    end
    
    % maybe do something 
    
%     pause(0.1)
end