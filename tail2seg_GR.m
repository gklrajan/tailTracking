%% tail2seg
% segment tail into 11 bits - use the first (last in the actual output) 9 bits for
% calculations for higher reliability.


%% Gokul Rajan | 2019-12-09 | to track tail in embeddded movies

% 2019/12/10 - modified to only largest blob detection
% 2019/12/10 - vis extracted coords on gaussian filt image
% 2019/12/11 - added a section to find angles of the 11 tail segs

% 2021/06/01 - checked and annotated with minor updates for Elena
% 2021/06/04 - ran a test - works well!

% tbd:
% time different sections

%%

clear; clc;
addpath(genpath('/Institut Curie/Lab/Projects/allScripts/ZebranalysisSystem/zebraFunctions/')); 
%rootDir='/Users/grajan/Desktop/';

%rootDir='/Institut Curie/Lab/Projects/Danionella/embeddedBehavior/ZF/converted/';
%rootDir='/Institut Curie/Lab/Projects/Danionella/embeddedBehavior/DT/6dpf+selectFinEmbedded+swim+20min_/converted/';
rootDir='/Institut Curie/Lab/Projects/Danionella/lightSheet/analysis/tail/DT/run01/';


%%
myDir = uigetdir('/Institut Curie/Lab/Projects/','Go to your directory!');
if ~isdir(myDir)
    uiwait(warndlg('The selected directory does not exist'));
    return;
end

filePattern = fullfile(myDir,'*.avi');
myFiles = dir(filePattern); %all txt files

cd(rootDir);

%%
for ff = 1:length(myFiles) % go through each file
    tic;
    
    dirName = myDir;
    fileName = myFiles(ff).name;
    inputName = fullfile(myDir,fileName);
    
    
v = VideoReader(fileName); % grab a video
kk=1;

vSample1 = read(v,70); % get 70th image
vSample2 = read(v,700); % get 700th image
vSample3 = read(v,7000); % get 7000th image

level1=thresh_tool(imadjust(rgb2gray(vSample1))); % GUI based thresh selection
level2=thresh_tool(imadjust(rgb2gray(vSample2)));
level3=thresh_tool(imadjust(rgb2gray(vSample3)));

level=round((level1+level2+level3)/3,0);
level=round((level/255),2); % averaged thresh value

repImage=imadjust(rgb2gray(vSample3));
repImage=imgaussfilt((repImage),2);
imshow(repImage); % look at imadjusted and gassian filt img

repImage= im2bw(repImage,level); % now apply the thresh
repImage=imcomplement(repImage); % get the complement img
[JJ, rect] = imcrop(repImage);
imshow(JJ); % check processed frame
disp('example frame used for processing');

duration=v.Duration;
numFrames=duration*v.FrameRate;

%% preallocate for speed
tailSegCoords_master=cell(numFrames,11);
tailSegAng_master=nan(numFrames,10);

%%
cc = VideoReader(fileName); % get a video
while hasFrame(cc)  % now go through each frame of the video one-by-one and process
        
    % process each frame
    frame = readFrame(cc);
    oriFrame=frame; % copy the raw frame
    frame=imadjust(rgb2gray(frame));
    frame= imcrop(frame,rect);
    frame=imgaussfilt((frame),2);
    oriFrameCropped=frame; % copy cropped adjusted original frame
    frame= im2bw(frame,level);
    frame=imcomplement(frame);
    
    % block-wise processing - this is a bummer as the output can't be a
    % cell/struct
    
%     [p,q]=size(frame);
%     m=p;
%     n=round((q-1)/12);
%     fun = @(block_struct) regionprops(block_struct.data,'centroid');    
%     blockCentroids(ii) = blockproc(frame,[m n],fun);
%     
%     clf;
%     imshow(frame); hold on;
%     plot(blockCentroids{ii},'*');
%     pause(1); hold off;
        
    
   %scan by region - reinventing a simple block scan to go over the image in
   %vertical sections
   % first apply imclose to join any nearby components in bin img
   se = strel('disk',5);
   frameClosed=  imclose(frame,se);
   
    [p,q]=size(frameClosed);
    m=p;
    n=round((q-1)/12); % divide the image into 12 vertical blocks
    
    %init
    frameRegs=cell(1,11);
    
    frameProcess=frameClosed; % now copy the frame to be processed
    mask = zeros(size(frameClosed)); % make a generic black mask
    mask(:,[1:n]) = 1; % adapt the mask for 1st block
        
    frameProcess = frameProcess.*mask; % img with only first block's binarized data!
    frameProcess= bwareafilt(logical(frameProcess),1); % select only the largest blob
    frameRegs{1}=regionprops(frameProcess,'centroid'); % find the centroid
    
    clearvars mask frameProcess;
    
    for ii=2:11 % now repeat the above process in the vertical blocks 2 to 11 - 12th seg near swim bladder is ignored
    
    frameProcess=frameClosed;        
    mask = zeros(size(frameClosed));
    mask(:,[n*(ii-1):n*(ii)]) = 1;
    frameProcess = frameProcess.*mask;
    frameProcess= bwareafilt(logical(frameProcess),1);
    frameRegs{ii}= regionprops(frameProcess,'centroid');
    %imshow(frameProcess); 
    clearvars mask frameProcess;
    
    end
    
     %frameRegsMat=nan(11,2);
    % storing the coords in a matrix
    for hh=1:size(frameRegs,2)
        if ~isempty(frameRegs{hh})
        frameRegsMat(hh,:)=frameRegs{hh}.Centroid;
        end
    end

    
    if ismember('frameRegsMat',evalin('base','who'))
        
    frameRegsMat=flip(frameRegsMat); % flip tail direction => now base to tip; use this to store final coords and angs
    
    tailSegAng=nan(1,10);
    % calculate angles of the 10 tail segments formed by 11 points
    for zz=1:(size(frameRegsMat,1)-1)
    tailSegAng(zz) = atan2(abs(det([frameRegsMat(zz,:);frameRegsMat(zz+1,:)])),dot(frameRegsMat(zz,:),frameRegsMat(zz+1,:)));
    end
        
    else
        frameRegsMat=nan(11,2);
        tailSegAng=nan(1,10);
    end % cal angle only if tail is segmented

    
    %%
        
%     %% visualize - show tracked points on the tail
%     imshow((oriFrameCropped)); hold on;
%     scatter(frameRegsMat(:,1),frameRegsMat(:,2),40,'MarkerEdgeColor',[0 .5 .5],...
%               'MarkerFaceColor',[0 .7 .7],...
%               'LineWidth',1.5); hold off;
%           %pause(0.2);

%%
    
   
    tailSegCoords_master{kk}=frameRegsMat;
    tailSegAng_master(kk,:)=tailSegAng;
    
    
    kk=kk+1;
    
    if rem(kk,10000)==0
        disp(kk);
    end
    
    clearvars -except cc kk tailSegCoords_master tailSegAng_master level rect fileName myDir myFiles filePattern rootDir duration;
    
end % end of processing a frame

toc; % time for processing one video

save(strcat(rootDir,fileName,'.mat'),'tailSegCoords_master','tailSegAng_master','duration','-v7.3');

clearvars -except myDir myFiles filePattern rootDir

end % end of processing a video
