%% tail2seg_process
% to process the extracted tail seg angles extracted using tail2seg

%% Gokul Rajan | 2019-12-11 | to process the extracted tail seg angles

% 2021-06-01 - checked for elena
% git versio ctrl started

%%
thresh=0.01; %0.25 for zf %0.02 for dt
numsegs=10;

%%

tailSegAngles=tailSegAng_master;
cumsumtailangles=cumsum(tailSegAngles')';

smoothedCumsumFixedSegmentAngles=cumsumtailangles;

%for the bout detection we smooth the tail curvature to eliminate kinks due
%to tracking noise
for n=2:size(cumsumtailangles,2)-1
    smoothedCumsumFixedSegmentAngles(:,n)=mean(cumsumtailangles(:,n-1:n+1)');
end

%we consider the difference in segment angles, because we want to detect
%tail movement
fixedSegmentAngles = [zeros(size(cumsumtailangles,2),1)'; diff(smoothedCumsumFixedSegmentAngles)];

%just take the number of well tracked segments
realSegmentAngles = fixedSegmentAngles(:,1:numsegs);
disp(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Filter angles%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this parameters set the timescale at which to smooth the tail movement
bcsize=10;
filteredSegmentAngle = realSegmentAngles*0;
for n = 1 : (size(realSegmentAngles,2))
    %smooth the tail movements at the shorter timescale
filteredSegmentAngle(:,n) = conv(realSegmentAngles(:,n),ones(bcsize,1)'/bcsize,'same');
end
disp(2)

%%%%%%%%%%%%%%%%%%%%%%%%%make cumsum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%we sum the angle differences down the length of the tail to give
%prominence to regions of continuous curvature in one direction
cumFilteredSegmentAngle = cumsum(filteredSegmentAngle')';

%%
%%%%%%%%%%%%%%%%%%%%%%%%%Calculate cumsum of cumsum%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%we sum up the absolute value of this accumulated curvature so bends in
%both directions are considered
superCumSegAngle = cumsum(abs(cumFilteredSegmentAngle)')';

%%
%%%%%%%%%%%%%%%%%%%%%%%%%Calculate tail curvature%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%it convolves all the segments in one
%this measure is filtered with a boxcar filter at the tail-beat timescale
bcFilt=10;
tailCurveMeasure = conv(superCumSegAngle(:,end),boxcarf(bcFilt),'same');
%%
%%%%%%%%%%%%%%%%%%%%%%min max filter to smooth data%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%min filter removes baseline fluctuations (needs to be well above the bout
%length)
%max filter flattens out beat variations (timescale needs to be adjusted
%according to be well below the interbout length)
smootherTailCurveMeasure = maxFilterFast(tailCurveMeasure(:,end),20) - minFilterFast(tailCurveMeasure(:,end),400);
%%
%%%%%%%%%%%%%%%%%%% fix size bug of smootherTailCurveMeasure %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(smootherTailCurveMeasure) ~= length(cumsumtailangles)
    smootherTailCurveMeasure(end) = [];
end


%%

% filter used on the curvature measure
filterB = ones(1,100)./100;

idx_nan3=isnan(smootherTailCurveMeasure);
smootherTailCurveMeasure_a=smootherTailCurveMeasure;

smootherTailCurveMeasure_a(idx_nan3) = 0;
smootherTailCurveMeasure_a_filtfilt        = filtfilt(filterB, 1, smootherTailCurveMeasure_a);
smootherTailCurveMeasure_a_filtfilt(idx_nan3) = nan; % re-insert the nan values
smootherTailCurveMeasure_a_filtfilt_intrp=interp1gap(smootherTailCurveMeasure_a_filtfilt,80,'linear');

%plot(smootherTailCurveMeasure); hold on; plot(smootherTailCurveMeasure_a_filtfilt_intrp);


%%
%I picked this threshold manually, but a wide range should work, so long as
%the data is not too noisy
allbout= smootherTailCurveMeasure_a_filtfilt_intrp>thresh; %0.25 for zf %0.02 for dt
idx_nan4=find(isnan(smootherTailCurveMeasure_a_filtfilt_intrp));

allstart=max([0 diff(allbout)],0);
allend=max([0 -diff(allbout)],0);

%forTBF=sgolayfilt(cumsumtailangles(:,numsegs),3,35); %8th seg - 2.4mm
%ff=rad2deg(cumFilteredSegmentAngle);

binVec20=smootherTailCurveMeasure;
binVec20(allbout)=1;
binVec20(~allbout)=0;


imagesc(binVec20);
% plot(cumsumtailangles(:,numsegs)); hold on; plot(forTBF);

%%
% ff=rad2deg(cumFilteredSegmentAngle);
% Fs=4;
% t = (0:length(ff)-1)/Fs;
% ff_tail=ff(:,8);
% 
% nanx = isnan(ff_tail);
% t    = 1:numel(ff_tail);
% int_ff_tail(nanx) = interp1(t(~nanx), ff_tail(~nanx), t(nanx));
% 
% h=imagesc(abs(ff_tail'));

%% save 10min bin vector - 10/10 mins
1;

binVec10=binVec20(1:150000);
imagesc(binVec10);


%% percent active time - 10/20 mins
% 10 min - 150000 - 150656
% 20 min - 300000 - 300562

%auto
%actTime_perCent=(sum(binVec20)./(length(binVec20)-length(idx_nan4))).*100;
%manual
actTime_perCent=(sum(binVec20)./(length(binVec20))).*100;

%% max duration - 10/20 mins

BinVec3D=[binVec20;binVec20;binVec20];
CC = bwconncomp(BinVec3D,4);
numPixels = cellfun(@numel,CC.PixelIdxList);
numPixels=numPixels+2; % for the trimmed 1 frame at onset and offset each - 'cos bwconncomp woeks only on 2d/3d


allSec=(numPixels.*1.43)./1000;

[biggest,idx] = nanmax(numPixels);
maxSec=(biggest.*1.43)./1000;

medSec = nanmedian(numPixels);
medSec=(medSec.*1.43)./1000;


clearvars -except allSec maxSec actTime_perCent binVec10 medSec binVec20 mast_allSec mast_maxSec mast_actTime_perCent mast_binVec10 mast_medSec mast_binVec20

% %%%%%%%%%%
% mast_binVec10=[];
% mast_allSec=[];
% mast_maxSec=[];
% mast_medSec=[];
% mast_actTime_perCent=[];
% %mast_binVec20=[];


%%%%%%%%%%%%
dcm1 = datacursormode(gcf);
set(dcm1, 'UpdateFcn', @Data_Cursor_precision, 'Enable', 'on');
1;
%%%%%%%%%%%%
mast_binVec10=[mast_binVec10;binVec10];
mast_allSec=[mast_allSec;allSec'];
mast_maxSec=[mast_maxSec;maxSec];
mast_medSec=[mast_medSec;medSec];
%mast_binVec20=[mast_binVec20;binVec20];
mast_actTime_perCent=[mast_actTime_perCent;actTime_perCent];

clearvars -except mast_allSec mast_maxSec mast_actTime_perCent mast_binVec10 mast_medSec mast_binVec20

1;
%%

%% plot scaled image
con=[mast_binVec10_DT;mast_binVec10_ZF];

con(con==1)=NaN;
con(con==0)=inf;

con(isnan(con))=0;
con(con==inf)=1;

imagesc(con(5:8,:));