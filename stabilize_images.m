% STABILIZEIMAGES()
% Stabilize images by removing jitter computed using phase correlation
% NPM 2019-2020
% 
% We compute the shifts to stabilize all volumes and write the 16
% bit stabilized volumes to disk. We generate the stabilized MIPs and show
% an RGB overlay of the reference timepoint data and stabilized data in
% cyan/magenta.
% 
%
% Requirements
% ------------
% data specified by fileName 
% This is accomplished by running make_mips.m
%
% Parameters
% ----------
% mipDir : where MIPs are stored
% mips_stab_check : output dir for the stabilized RGB overlays 
% mipoutdir : output dir for the stabilized MIPs
% im_intensity : scale for MIP output intensity
% imref_intensity : scale for reference MIP intensity in RGB overlay
% t_ref : timestamp of reference, matches other timepoints to this one
% timePoints : all timestamps to consider in the sequence
% times_todo : timestamps to correct/make output
% typename : output volume datatype, as string ('uint8', 'uint16')
% fileName : input filename format string for data to stabilize
% fileNameOut : output filename format string for stabilized volumes
% rgbName : RGB overlay filename format string
%
% Returns
% -------
% - jitter/drift-corrected 3D volumes 
% - MIPs of the jitter/drift-corrected volumes within mipoutdir
% - colored overlays of one view with the reference time. 

%% Make the subdirectories for the mips if not already existing
t_ref_ind = find( timePoints == t_ref ) ;

%% Define shifts for each time point ======================================
disp('Defining shifts...')


disp('Shifts not on disk, computing them...')
% Load MIP data into im_1 and im_2 for all times
disp('Loading MIP data for all times...')
NTimes = length(timePoints);
% preallocate im_1 for speed
tmp = imread(fullfile(mipDir, sprintf(name1,timePoints(1)))) ;
im_1 = zeros([size(tmp) length(timePoints)]) ;
% preallocate im_2 for speed
tmp = imread(fullfile(mipDir, sprintf(name11,timePoints(1)))) ;
im_2 = zeros([size(tmp) length(timePoints)]) ;
for tid = 1:length(timePoints)
    time = timePoints(tid) ;
    % weight the different views equally
    im1a = fullfile(mipDir, sprintf(name1, time)) ;
    im1b = fullfile(mipDir, sprintf(name2, time)) ;
    im_1(:,:,tid) = imread(im1a) + imread(im1b) ;

    im2a = fullfile(mipDir, sprintf(name11, time)) ;
    im2b = fullfile(mipDir, sprintf(name21, time)) ;
    im_2(:,:,tid) = imread(im2a) + imread(im2b) ;

    im3a = fullfile(mipDir, sprintf(name12, time)) ;
    im3b = fullfile(mipDir, sprintf(name22, time)) ;
    im_3(:,:,tid) = imread(im3a) + imread(im3b);
end
disp('done loading data into im_1 and im_2')

% Compute shifts via phase correlation
disp('Computing/overwriting shifts')
shifts = struct('x_1',[],'y_1',[],'x_2',[],'y_2',[]);
for time = 1 :NTimes
    [shiftx,shifty,~] = xcorr2fft(im_1(:,:,time), im_1(:,:,t_ref_ind));
    shifts(time).x_1 = shiftx;
    shifts(time).y_1 = shifty; 
end
disp('done defining phase correlations')

%% Convert correlations to shifts
dx =  cat(1,shifts.x_1); % rows        (x) 
dy =  cat(1,shifts.y_1); % columns     (y)

%% Plot the shifts
disp('Plotting shifts...')
close('all')
hold all;
% Consider view 0
plot(timePoints, dx, '.-', 'DisplayName', 'dx')
plot(timePoints, dy, '.-', 'DisplayName', 'dy')
ylabel('shift [pixels]')
xlabel('timestamp')
legend('location', 'best')
title('Jitter stabilization')
saveas(gcf, fullfile(mipDir, 'jitter_stabilization.png'))
disp('done plotting shifts, see Figure')
disp('Saving shifts to shifts_stab.mat in mipDir')
save(fullfile(mipDir, 'shifts_stab.mat'), 'shifts')


close('all')


%% Build reference MIP for RGB overlay
disp('building reference MIP...')
name_ref = sprintf(fileName, t_ref);
% preallocate im_ref3D for speed
im_ref3D = zeros([size(tmp) stackSize]) ;
for z = 1 : stackSize
    im_ref3D(:,:,z)= imread(name_ref,z);
end
mip_ref = squeeze(max(im_ref3D,[],3));
clear tmp
disp('done creating reference MIP')

%% Build image for each timepoint =========================================
disp('Running through timepoints to build ims...')
for tid = 1 : length(timePoints)
    disp(['considering tid = ' num2str(tid)])
    time = timePoints(tid);
    if ismember(time, times_todo)
        % The original image in 3d
        im0fn = sprintf(fileName, time);

        % Check that we're not appending to an existing file
        name_out = sprintf(fileNameOut,time) ;
        tiff_exists = exist(name_out, 'file') ;
        if tiff_exists && ~overwrite_tiffs && ~overwrite_mips
            disp(['Output file already exists: ' name_out ])
        else
            disp('Creating mips for this timepoint')
            % preallocate im_2_3D for speed, and 
            % Specify typename for correct output bitdepth
            tmp = imread(im0fn, 1) ;
            im0 = zeros([size(tmp) stackSize], typename) ;
            clear tmp

            for z = 1:stackSize
                im0(:,:,z)= imread(im0fn,z);
            end

            % Initialize the image 
            im = im0 ;
            imx = 0 * im0 ;  % for shift in x

            % Offset if there is a shift in X
            if dx(tid)~=0
                if dx(tid)>0
                    % imx(dx(tid):end,:,:) = im(1:end-dx(tid)+1,:,:);
                    imx(1+dx(tid):end,:,:) = im(1:end-dx(tid),:,:);
                else
                    % imx(1:(end+dx(tid)+1),:,:) = im(-dx(tid):end,:,:);
                    imx(1:(end+dx(tid)),:,:) = im((1-dx(tid)):end,:,:);
                end
            else
                imx = im ;
            end

            % Offset if there is a shift in Y
            if dy(tid)~=0
                if dy(tid)>0
                    % imy(:,dy(tid):end,:) = imx(:,1:(end-dy(tid)+1),:);
                    im(:,1+dy(tid):end,:) = imx(:,1:(end-dy(tid)),:);
                else
                    % imy(:,1:(end+dy(tid)+1),:) = imx(:,-dy(tid):end,:);
                    im(:,1:(end+dy(tid)),:) = imx(:,(1-dy(tid)):end,:);
                end
            else
                im = imx ;
            end
            
            % Make a color of the current mip wrt the reference
            mip_1 = squeeze(max(im,[],3));
            if strcmp(typename, 'uint8')
                rgb = zeros(size(mip_ref,1),size(mip_ref,2),3,'uint8');
                rgb(:,:,1)= uint8(mip_1 * im_intensity);
                rgb(:,:,2)= uint8(mip_1 * im_intensity);
                rgb(:,:,2)= uint8(mip_ref * im_intensity);
                rgb(:,:,3)= uint8(mip_ref * imref_intensity);
                rgb(rgb(:) > 255) = 255 ;
            elseif strcmp(typename, 'uint16')
                rgb = zeros(size(mip_ref,1),size(mip_ref,2),3,'uint16');
                rgb(:,:,1)= uint16(mip_1 * im_intensity);
                rgb(:,:,2)= uint16(mip_1 * im_intensity);
                rgb(:,:,2)= uint16(mip_ref * im_intensity);
                rgb(:,:,3)= uint16(mip_ref * imref_intensity);
                rgb(rgb(:) > 65535) = 65535 ;
            else
                error(['Have not coded for RGB overlay ', ...
                    'with given typename: ', typename, '. Do so here'])
            end
            
            % Make record of the mip drift
            imwrite(rgb, fullfile(mipsRGBDir, sprintf(rgbName,time)))
        end
    end
end
disp('done building ims and saving stab tifs...')

% 
