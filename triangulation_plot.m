load ('position_list_1.mat')
filenamestr = "../20191017181230_JGFPHRFP_1um_10spf_2x63x_oiloil_4pc1pc_max_z_RFP.tif";
% 0.361984 micron/pixel

m_length = zeros(1,size(objects_positions,1)-1);
std_length = m_length;
std_err_length = m_length;


for i =  1: size(objects_positions,1)-1
    
    save_filename = strcat('triangulation/triangulation_', sprintf('%03d',i), '.png');
    save_filename_tris = strcat('triangles/triangulation_', sprintf('%03d',i), '.mat');
    img1 = imread(filenamestr, i);
    P = objects_positions{i,1}.*size(img1); % points
    %px = P(:, 1);
    %py = P(:, 2);
    DT = delaunayTriangulation(P);  
    %ratio_list = metriclist(DT.ConnectivityList, P);
    cl = clean_connlist(0.87, DT.ConnectivityList, P);
    TR = triangulation(cl, P);
    edge_list = edges(TR);
    length_list = get_length(edge_list, P);
    length_list = length_list*0.361984;
    F = freeBoundary(TR);
    
    m_length(i) = mean(length_list);
    std_length(i) = std(length_list);
    std_err_length(i) = std_length(i)/sqrt(size(length_list,1));
    subplot(1,2,1)
    triplot(TR);
    hold on
    
%    % plot txt
%     for k = 1: size(DT.ConnectivityList,1)
%         text(ratio_list(k, 2), ratio_list(k, 3), num2str(ratio_list(k, 1)), 'Color', 'r')
%     end
    
%   %plot boundary
%   plot(px(F),py(F),'-r','LineWidth',2);

    imh = imshow(img1);
    hold off
    uistack(imh, 'bottom')

    subplot(1,2,2)
    histogram(length_list,20)
    xlim([0 50])
    xlabel("distance [micron]")
    ylim([0 150])
    %saveas(gcf,save_filename)
    
    save(save_filename_tris, 'TR')
end




close(gcf)
%%

save('lengths.mat','m_length', 'std_length', 'std_err_length')
time = 0:10:(length(m_length)-1)*10;
t2 = [time, fliplr(time)];

hold on;
std_plus = m_length+std_length;
std_minus = m_length-std_length;

std_err_plus = m_length+std_err_length;
std_err_minus = m_length-std_err_length;

inBetween1 = [std_minus, fliplr(std_plus)];
inBetween2 = [std_err_minus, fliplr(std_err_plus)];
h1 = fill(t2, inBetween1, 'r');
h2 = fill(t2, inBetween2, 'g');
set(h1,'facealpha',.5)
set(h2, 'facealpha', .5)
plot(time ,m_length, 'Color', 'b', 'LineWidth', 2);
xlabel("time [s]")
ylabel("average distance [micron]")

%%

function length_list = get_length(edge_list, point_list)
    length_list = zeros(size(edge_list,1),1);
    for i = 1: size(edge_list, 1)
        pos1 = point_list(edge_list(i,1), 1);
        pos2 = point_list(edge_list(i,2), 1);
        distance = sqrt(dot(pos1-pos2, pos1-pos2));
        length_list(i) = distance;
    end
end


function cl = clean_connlist(cutoff, connlist, point_list)
    for i = 1: size(connlist, 1)
        p1 = point_list(connlist(i, 1), :);
        p2 = point_list(connlist(i, 2), :);
        p3 = point_list(connlist(i, 3), :);
        ratio = reject_triangle(p1,p2,p3);
        if ratio > cutoff
            connlist(i, :) = NaN(1,3);
        end
    end
    cl = rmmissing(connlist);
end




function ratio_list = metriclist(connlist, point_list)
    ratio_list = zeros(size(connlist, 1), 3);
    for i = 1: size(connlist, 1)
        p1 = point_list(connlist(i, 1), :);
        p2 = point_list(connlist(i, 2), :);
        p3 = point_list(connlist(i, 3), :);
        center = (p1+p2+p3)/3;
        ratio_list(i,1) = reject_triangle(p1,p2,p3);
        ratio_list(i,2:3) = center;
    end
end
