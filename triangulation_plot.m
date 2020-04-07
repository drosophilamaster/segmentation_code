load ('position_list_1.mat')
filenamestr = "../20191017181230_JGFPHRFP_1um_10spf_2x63x_oiloil_4pc1pc_max_z_RFP.tif";


m_length = zeros(1,size(objects_positions,1)-1);
for i = 1:size(objects_positions,1)-1
    
    save_filename = strcat('triangulation/triangulation_', sprintf('%03d',i), '.png');
    img1 = imread(filenamestr, i);
    P = objects_positions{i,1}.*size(img1);
    px = P(:, 1);
    py = P(:, 2);
    DT = delaunayTriangulation(P);
    DT.ConnectivityList
    TR = triangulation(DT.ConnectivityList, P);
    edge_list = edges(TR);
    length_list = get_length(edge_list, P);
    F = freeBoundary(TR);
    
    m_length(i) = mean(length_list);
    
    
    
    subplot(1,2,1)
    triplot(TR);
    hold on
    plot(px(F),py(F),'-r','LineWidth',2);
    imh = imshow(img1);
    hold off
    uistack(imh, 'bottom')

    subplot(1,2,2)
    histogram(length_list,20)
    xlim([0 200])
    ylim([0 250])
    saveas(gcf,save_filename)
end

plot([0:10:(length(m_length)-1)*10],m_length);
xlabel("time [s]")
ylabel("average distance [pix]")



function length_list = get_length(edge_list, point_list)
    length_list = zeros(size(edge_list,1),1);
    for i = 1: size(edge_list, 1)
        pos1 = point_list(edge_list(i,1), 1);
        pos2 = point_list(edge_list(i,2), 1);
        distance = sqrt(dot(pos1-pos2, pos1-pos2));
        length_list(i) = distance;
    end
end