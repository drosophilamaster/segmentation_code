filenamestr = "../20191017181230_JGFPHRFP_1um_10spf_2x63x_oiloil_4pc1pc_max_z_RFP.tif";
num_timepoints = size(imfinfo(filenamestr),1);
pos_x = 75; pos_y = 75;
for i = 1:num_timepoints-2
    big_img1 = zeros(600,400);
    big_img2 = zeros(600, 400);
    img1 = imread(filenamestr, i);
    img2 = imread(filenamestr, i+1);
    [shiftx, shifty, c] = xcorr2fft(img1, img2);
    big_img1(pos_x+shiftx:pos_x+shiftx+255,pos_y+shifty:pos_y+shifty+255) = img1;
    big_img2(pos_x:pos_x+255, pos_y:pos_y+255) = img2;
    pos_x = pos_x-shiftx;
    pos_y = pos_y-shifty;
    imshowpair(big_img1, big_img2);
end



