clear;

load ('position_list_1.mat')
filenamestr = "../20191017181230_JGFPHRFP_1um_10spf_2x63x_oiloil_4pc1pc_max_z_RFP.tif";

mean_psi6 = zeros(1,size(objects_positions,1)-1);
mean_mod_psi = zeros(1,size(objects_positions,1)-1);

for i = 1:size(objects_positions,1)-1
    
    save_filename = strcat('psi6_plots/psi6_', sprintf('%03d',i), '.png');
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
    
    
    
    psi1 = rmmissing(psi);
    mean_mod_psi(i) = mean(abs(psi1));
    
% 
%     
%     subplot(1,2,1)
%     triplot(TR);
%     hold on
%     
%     % plot txt
%      for k = 1: num_points
%          text(points(k, 1), points(k, 2), num2str(psi(k)), 'Color', 'r')
%      end
%     
% 
%     imh = imshow(img1);
%     hold off
%     uistack(imh, 'bottom')
% 
%     subplot(1,2,2)
%     scatter(real(psi1), imag(psi1))
%     xlim([-1.1, 1.1])
%     ylim([-1.1, 1.1])
%     axis equal    
%     xlim([-1.1, 1.1])
%     ylim([-1.1, 1.1])
%     hold on;
%     plot(cos(0:0.01:2*pi), sin(0:0.01:2*pi), 'k-')
%     xlabel(['\Re \psi_' num2str(6)])
%     ylabel(['\Im \psi_' num2str(6)])
%     %saveas(gcf,save_filename)
%     hold off
end






function angle = get_angle(p0, p1)
    u = p1-p0;
    angle = atan(u(2)/u(1));
end