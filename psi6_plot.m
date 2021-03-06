clear;

load ('position_list_1.mat')
filenamestr = "../20191017181230_JGFPHRFP_1um_10spf_2x63x_oiloil_4pc1pc_max_z_RFP.tif";
load('psi6.mat')
cmap = colormap;
% mean_psi6 = zeros(1,size(objects_positions,1)-1);
% m_psi = zeros(1,size(objects_positions,1)-1);
% std_psi = m_psi;
% stderr_psi = m_psi;
for i = 1%1:size(objects_positions,1)-1
    
    save_filename = strcat('psi6plots/psi6_', sprintf('%03d',i), '.png');
    img1 = imread(filenamestr, i);
   
    load(strcat('triangles/triangulation_', sprintf('%03d',i), '.mat'));
    boundary_edges = TR.freeBoundary;
    boundary_points = boundary_edges(:, 1); 
    points = TR.Points;
    
    num_points = size(points, 1);
    conn_matrix = zeros(size(points, 1));
    
    for k = 1: num_points
        for kk = 1: k-1
            if isConnected(TR, k, kk)
                conn_matrix(k,kk) = 1;
            end
        end
    end
    
    
    conn_matrix = conn_matrix' + conn_matrix;

    psi = zeros(num_points, 1);
    
    
    for k1 = 1: num_points
        if ismember(k1, boundary_points)
            psi(k1) = nan();
        else
            points_attached = find(conn_matrix(k1, :) == 1);
            sum_exponential = 0;
            for k2 = 1:size(points_attached, 2)
                points(k1, :);
                points(points_attached(k2), :);
                angle = get_angle(points(k1, :), points(points_attached(k2), :));
                sum_exponential = sum_exponential + exp(6*1j*angle);
            end
            psi(k1) = sum_exponential/size(points_attached, 2);
        end
    end     
    
    
%     %% get mean, std, std err
     psi1 = rmmissing(psi);
%     m_psi(i) = mean(abs(psi1));
%     std_psi(i) = std(abs(psi1));
%     stderr_psi(i) = std(abs(psi1))/sqrt(length(psi1));
%     
    grayscaleIM(:,:,3) = img1;
    grayscaleIM(:,:,2) = img1;
    grayscaleIM(:,:,1) = img1;

    clf;
    subplot(1,2,1)
    imagesc(grayscaleIM);
    axis equal
    axis off
    hold on
   
    %hold off
    triplot(TR);
    hold on
    xs = points(:,1);
    ys = points(:,2);
    [v, c] = voronoin(points);
    %voronoi(xs,ys)
%     
%     hold on
   
    % plot txt
%      for k = 1: num_points
%          text(points(k, 1), points(k, 2), num2str(psi(k)), 'Color', 'r')
%      end
%     
    for k = 1: num_points
        psi_value = psi(k);
        if isnan(psi_value)
        else
            fill_color(abs(psi_value), cell2mat(c(k,1)), v, cmap)
        end
        
    end
    

    % Add colorbar
    c = colorbar() ;

    hold on
    %imh = imshow(img1); 
    %uistack(imh, 'bottom')
    

    
    subplot(1,2,2)
    time = 0:10:(length(m_psi)-1)*10;
    t2 = [time, fliplr(time)];

    hold on;
    std_plus = m_psi+std_psi;
    std_minus = m_psi-std_psi;

    std_err_plus = m_psi+stderr_psi;
    std_err_minus = m_psi-stderr_psi;

    inBetween1 = [std_minus, fliplr(std_plus)];
    inBetween2 = [std_err_minus, fliplr(std_err_plus)];
    h1 = fill(t2, inBetween1, 'y');
    h2 = fill(t2, inBetween2, 'g');
    set(h1,'facealpha',.5)
    set(h2, 'facealpha', .5)
    plot(time ,m_psi, 'Color', 'b', 'LineWidth', 2);
    plot(time(i), m_psi(i), 'rs')
    ylim([0.1,0.8])
    xlim([0,3500])
    
    xlabel("time [s]")
    ylabel(['\psi_' num2str(6)])
    
    %set(gcf, 'Position', [100,100,1500,1000])
    
%     subplot(1,2,3)
%     hold on
%     cinds = uint8(min(abs(psi1),1)*length(cmap));
%     scatter(real(psi1), imag(psi1), 40, cmap(max(1, cinds),:));
%     xlim([-1.1, 1.1])
%     ylim([-1.1, 1.1])
%     axis equal    
%     xlim([-1.1, 1.1])
%     ylim([-1.1, 1.1])
%     plot(cos(0:0.01:2*pi), sin(0:0.01:2*pi), 'k-')
%     xlabel(['\Re \psi_' num2str(6)])
%     ylabel(['\Im \psi_' num2str(6)])
    %set(gcf, 'PaperUnits', 'centimeters');
    %set(gcf, 'PaperSize', [100 16]);
    set(gcf, 'Position', [100,100,1000,400])
    saveas(gcf,save_filename)
    hold off
end

%%
save('psi6.mat','m_psi', 'std_psi', 'stderr_psi')

%%
% load('lengths.mat')
% for i = flip(1:20, 2)%1:size(objects_positions,1)-1
%     save_filename = strcat('psi1/psi6_', sprintf('%03d',i), '.png');
%     load(strcat('triangles/triangulation_', sprintf('%03d',i), '.mat'))
%     img1 = imread(filenamestr, i);
%     clf
%     subplot(2,2,2)
%     triplot(TR)
%     
%     hold on;
%     
%     imh = imshow(img1);
%     hold off
%     uistack(imh, 'bottom')
%     
%     
%     subplot(2,2,1)
%     %% plot mean
%     time = 0:10:(length(m_psi)-1)*10;
%     t2 = [time, fliplr(time)];
% 
%     hold on;
%     std_plus = m_psi+std_psi;
%     std_minus = m_psi-std_psi;
% 
%     std_err_plus = m_psi+stderr_psi;
%     std_err_minus = m_psi-stderr_psi;
% 
%     inBetween1 = [std_minus, fliplr(std_plus)];
%     inBetween2 = [std_err_minus, fliplr(std_err_plus)];
%     h1 = fill(t2, inBetween1, 'y');
%     h2 = fill(t2, inBetween2, 'g');
%     set(h1,'facealpha',.5)
%     set(h2, 'facealpha', .5)
%     plot(time ,m_psi, 'Color', 'b', 'LineWidth', 2);
%     plot(time(i), m_psi(i), 'rs')
%     ylim([0.1,0.8])
%     xlim([0,3500])
% %     axis equal
% %     ylim([0.1,0.8])
% %     xlim([0,3500])
%     
%     xlabel("time [s]")
%     ylabel(['\psi_' num2str(6)])
%     
%     subplot(2,2,3)
%     hold on;
%     std_plus = m_length+std_length;
%     std_minus = m_length-std_length;
%     
%     std_err_plus = m_length+std_err_length;
%     std_err_minus = m_length-std_err_length;
%     
%     inBetween1 = [std_minus, fliplr(std_plus)];
%     inBetween2 = [std_err_minus, fliplr(std_err_plus)];
%     h1 = fill(t2, inBetween1, 'y');
%     h2 = fill(t2, inBetween2, 'g');
%     set(h1,'facealpha',.5)
%     set(h2, 'facealpha', .5)
%     plot(time ,m_length, 'Color', 'b', 'LineWidth', 2);
%     plot(time(i), m_length(i), 'rs')
%     xlabel("time [s]")
%     ylabel("average distance [micron]")
%     
%    
%     set(gcf, 'Position', get(0, 'Screensize'))
%     %pause(2)
%     saveas(gcf,save_filename)
%     hold off;
%end






function fill_color(value, indices, vertices, cmap)
    x = vertices(indices, 1);
    y = vertices(indices, 2);
    cind = uint8(min(1, value)*length(cmap));
    if max(x)<=256 && max(y) <= 256 && min(x)>=0 && min(y) >=0
        fill(x,y,  cmap(max(cind, 1), :), 'FaceAlpha', 0.5)
    end
end
function angle = get_angle(p0, p1)
    u = p1-p0;
    angle = atan(u(2)/u(1));
end