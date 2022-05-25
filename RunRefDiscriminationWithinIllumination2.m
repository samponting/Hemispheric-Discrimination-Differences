clear 
close all

shape = 'Bumpy';
specular = 'Matte';
environment = 'E1';
observer = 'Test';
session = 1;

% number of reference angle
noa = 8;

% number of reflectance
nor = 30;

% number of practice session
nop =  2;

conditionname = strcat(shape,'_',specular,'_',environment,'_');

addpath(genpath('C:\Data\Experiments\Reflectance Discrimination'));
rmpath('C:\Data\Experiments\Reflectance Discrimination\RGB_FixedView');
%rmpath('C:\Data\Experiments\Reflectance Discrimination\Stimulus RGB');

%cd('C:\Data\Experiments\Reflectance Discrimination');

load('ControlCondition_Chromaticity');
load('MB_AllHue.mat');
%load('noise_1pixel.mat');
%load('noise_4pixels.mat');
load('noise16pixels.mat');
%noise = noise_16pixel;
%number of hues tested
noh = 8;
stop = zeros(noh,1);
stop2 = zeros(noh,1);
%%%%%%%%%%%%%%%%%%%  Palameter Setting For Palamedes  %%%%%%%%%%%%%%%%%%%  
%Set up psi
NumTrials = 40;

grain = 101; %grain of posterior, high numbers make method more precise at the cost of RAM and time to compute.
             %Always check posterior after method completes [using e.g., :
             %image(PAL_Scale0to1(PM.pdf)*64)] to check whether appropriate
             %grain and parameter ranges were used.
             
PF = @PAL_Gumbel; %assumed psychometric function

%Stimulus values the method can select from
stimRange = 1:nor;

%Define parameter ranges to be inclPMed in posterior
priorAlphaRange = linspace(1,nor,grain);
priorBetaRange =  linspace(log10(.0625),log10(10),grain); %Use log10 transformed values of beta (slope) parameter in PF
priorGammaRange = 0.25;  %fixed value (using vector here would make it a free parameter) 
priorLambdaRange = .02; %ditto

%Initialize PM structure
for i = 1:noh
    PM(i) = PAL_AMPM_setupPM('priorAlphaRange',priorAlphaRange,...
                          'priorBetaRange',priorBetaRange,...
                          'priorGammaRange',priorGammaRange,...
                          'priorLambdaRange',priorLambdaRange,...
                          'numtrials',NumTrials,...
                          'PF' , PF,...
                          'stimRange',stimRange); 
    PM(i).xCurrent = nor;
end
PM2 = PM;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stimulus duration [s]
StimulusDuration = 0.2;
AdaptationDuration = 10;

% Blank [s]
PauseAfterTrial = 0;

%trial counter for each hue
Trial_cnt = ones(noh,1);
Trial_cnt2 = ones(noh,1);

endflag = 0;

%%%%%%%%%%%%%%%% Initilize VSG
disp('Initialising VSG...')
global CRS
crsInit('C:\Data\Equipment\ViSaGe + NEC CRT\Monitor calibration\190617\19June2017_1600x1200.vsg');
framerate = crsGetFrameRate;
crsSetVideoMode(CRS.HYPERCOLOURMODE); %42bit

crsResponseBoxOpen(CRS.respCEDRUS);
disp('...done');

blankPage = crsGetSystemAttribute(CRS.NUMVIDEOPAGES);
crsSetDisplayPage(blankPage);
screenWidth = crsGetScreenWidthPixels;
screenHeight = crsGetScreenHeightPixels;

ViewingDistance = 92; %cm

% pixel per degree
ppd = round(screenWidth/(2*atand((43/2)/ViewingDistance)));

%gap between images
%gapsize = round(0.49*ppd); %degree * ppd
gapsize = round(0.75*ppd);

screenPos=[0 0];
coordsTarget = screenPos(1,:);

for i=1:crsGetSystemAttribute(CRS.NUMVIDEOPAGES)
    crsSetDrawPage(CRS.VIDEOPAGE,i,[0 0 0]);
end
%sounds
t0 = zeros(1,0.1*44100);
ts=1:(44100*0.2);
s15=[t0 0.1*sin(ts/44100*2*pi*1500)];
sx=[t0 0.1*(sin(ts/44100*2*pi*800)+sin(ts/44100*2*pi*400)/2+sin(ts/44100*2*pi*200)/4+sin(ts/44100*2*pi*100)/8)];

disp('Beginning experiment')
makeChromTemporalNoiseMask_TM(2,noise);

disp('Press centre button to begin');
pause(AdaptationDuration);
sound(s15,44100);

WaitForCentreButton;
crsSetCommand(CRS.CYCLEPAGEDISABLE);
crsSetZoneDisplayPage(CRS.VIDEOPAGE,1);
tic;

%HideCursor;

practice = 1;

%start trials loop
while ~(sum(stop) == noh && sum(stop2) == noh) 
for h = randperm(noh)
    for Hemifield = randperm(2)
        if Hemifield == 1
            if PM(h).stop == 0
            degree = (h-1)*360/noh;
            degree
            Angle_seed = 45*(randperm(noa)-1);
            RefAngle1 = Angle_seed(1);
        %     RefAngle2 = Angle_seed(2);
        %     RefAngle3 = Angle_seed(3);
        
        %     file_Ref1 =  strcat(conditionname,'Ref_angle',num2str(RefAngle1),'_RGB.mat');
        %     file_Ref2 =  strcat(conditionname,'Ref_angle',num2str(RefAngle2),'_RGB.mat');
        %     file_Ref3 =  strcat(conditionname,'Ref_angle',num2str(RefAngle3),'_RGB.mat');
            file_Ref1 =  strcat(pwd,'/Stimuli/Img_ang',num2str(RefAngle1),'_n1.mat');
        %     file_Ref2 =  strcat('/Stimuli/Img_ang',num2str(RefAngle2),'_n1.mat');
        %     file_Ref3 =  strcat('/Stimuli/Img_ang',num2str(RefAngle3),'_n1.mat');
            %file_Ref1 =  strcat(conditionname,num2str(degree),'deg_n1_RGB.mat');
            %file_Ref2 =  strcat(conditionname,num2str(degree),'deg_n1_RGB.mat');
            %file_Ref3 =  strcat(conditionname,num2str(degree),'deg_n1_RGB.mat');
            PM(h).xCurrent
            file_Test =  strcat(pwd,'/Stimuli/Img_ang',num2str(degree),'_n',num2str(PM(h).xCurrent),'.mat');
        
            %file_Ref =  strcat('UniformPatch_',num2str(degree),'deg_n1.mat');
            %file_Test =  strcat('UniformPatch_',num2str(degree),'deg_n',num2str(Intensity(h,Trial_cnt(h,1))),'.mat');
            
            factor = 3.2117*1.5/4.7;
            % 225.5 pixel 3 degree   
            % 188.9 pixel 2.5 degree   
        
            load(file_Test);Image(:,:,:,1) = stim;
            load(file_Ref1);Image(:,:,:,2) = stim;
        %     load(file_Ref2);Image(:,:,:,3) = stim;
        %     load(file_Ref3);Image(:,:,:,4) = stim;
        
            %file_TestName{Trial,h} = file_Test;
        
            imsize = size(Image);
        
            % Image to display stimuli.
            StimuliImage = zeros(imsize(2)*2+gapsize,imsize(1)*2+gapsize,3);
        
            arrangement = randperm(2);
            ansIndex(h,Trial_cnt(h,1)) = find(arrangement==1);
        
            % n
            % 1 3
            % 2 4
            r1 = 1:imsize(1);c1 = 1:imsize(2);
            r2 = imsize(1)+2*gapsize:2*imsize(1)+2*gapsize-1;c2 = 1:imsize(2);
            r3 = 1:imsize(1);c3 = imsize(2)+2*gapsize:2*imsize(2)+2*gapsize-1;
            r4 = imsize(1)+2*gapsize:2*imsize(1)+2*gapsize-1;c4 = imsize(1)+2*gapsize:2*imsize(1)+2*gapsize-1;
            
            StimuliImage(r1,c1,:) = Image(:,:,:,arrangement(1));
            StimuliImage(r3,c3,:) = Image(:,:,:,arrangement(2));
        %     StimuliImage(r3,c3,:) = Image(:,:,:,arrangement(3));
        %     StimuliImage(r4,c4,:) = Image(:,:,:,arrangement(4));
        
            StimuliImage = max(StimuliImage,0);
            StimuliImage = min(StimuliImage,1);
            Imsize = size(StimuliImage);
        
            crsSetDrawPage(CRS.VIDEOPAGE,1,1);
            crsSetCommand(CRS.CYCLEPAGEDISABLE);
            crsDrawMatrix42bitColour(coordsTarget,double(StimuliImage(:,1:2:Imsize(2),:)));
            crsSetZoneDisplayPage(CRS.VIDEOPAGE,1);
            
            pause(StimulusDuration);
            
            makeChromTemporalNoiseMask_TM(2,noise);
            crsResponseBoxFlush;
        
            button = -1;
            while button == -1;
                BR = crsGetResponse(CRS.respCEDRUS);
                pause(0.01);
        %         if (BR.Top == 1)||(BR.Left == 1)||(BR.Bottom == 1)||(BR.Right == 1)
        %             button = 1;
        %         end
                if (BR.Left == 1)||(BR.Right == 1)
                    button = 1;
                end
            pause(0.01);
            end
        
        %     if (BR.Top == 1) % top button
        %         response(h,Trial_cnt(h,1)) = 1;
        %     end
            if (BR.Left == 1) % left button
                response(h,Trial_cnt(h,1)) = 2;
            end
        %     if (BR.Bottom == 1) % bottom button
        %         response(h,Trial_cnt(h,1)) = 4;
        %     end
            if (BR.Right == 1) % right button
                response(h,Trial_cnt(h,1)) = 3;
            end
        
            %calculate performance or detect forced exit
            if response(h,Trial_cnt(h,1)) == ansIndex(h,Trial_cnt(h,1))
                sound(s15,44100);
                if practice == 0
                    correct(h,Trial_cnt(h,1)) = 1;
                end
            elseif response(h,Trial_cnt(h,1)) ~= ansIndex(h,Trial_cnt(h,1))
                sound(sx,44100);
                if practice == 0
                    correct(h,Trial_cnt(h,1)) = 0;
                end
            end
            pause(PauseAfterTrial);
        
            % Update intensity based on palamedes
            if practice == 0
                PM(h) = PAL_AMPM_updatePM(PM(h), correct(h,Trial_cnt(h,1))); %update PM structure
                % end if last 5 responses are stable
                if Trial_cnt(h,1) > 10 && mean(PM(h).seThreshold(Trial_cnt(h,1)-5:Trial_cnt(h,1)-1)) < 2
                    PM(h).stop = 1;
                end
            
                %if (PM(1).stop)&&(PM(2).stop)&&(PM(3).stop)&&(PM(4).stop)&&(PM(5).stop)&&(PM(6).stop)&&(PM(7).stop)&&(PM(8).stop) == 1
                
                for i = 1:noh
                    stop(i,1) = PM(i).stop;
                end
                
                if (sum(stop2) == noh && sum(stop) == noh)
                    pause(0.5);
                    [tada, Fs] = wavread('C:\Windows\Media\tada.wav');
                    wavplay(2*tada,Fs);
                    vsgBlack
                end
            end
            
            Trial_cnt(h,1) = Trial_cnt(h,1)+1;
            
            if min(Trial_cnt(:,1) == (nop+1)*ones(noh,1))==1 && practice == 1
                Trial_cnt(:,1) = ones(noh,1);
                practice = 0;
            end
            
            end
        elseif Hemifield == 2
            if PM2(h).stop == 0
                degree = (h-1)*360/noh;
                degree
                Angle_seed = 45*(randperm(noa)-1);
                RefAngle1 = Angle_seed(1);
            %     RefAngle2 = Angle_seed(2);
            %     RefAngle3 = Angle_seed(3);
            
            %     file_Ref1 =  strcat(conditionname,'Ref_angle',num2str(RefAngle1),'_RGB.mat');
            %     file_Ref2 =  strcat(conditionname,'Ref_angle',num2str(RefAngle2),'_RGB.mat');
            %     file_Ref3 =  strcat(conditionname,'Ref_angle',num2str(RefAngle3),'_RGB.mat');
                file_Ref1 =  strcat(pwd,'/Stimuli/Img_ang',num2str(RefAngle1),'_n1.mat');
            %     file_Ref2 =  strcat('/Stimuli/Img_ang',num2str(RefAngle2),'_n1.mat');
            %     file_Ref3 =  strcat('/Stimuli/Img_ang',num2str(RefAngle3),'_n1.mat');
                %file_Ref1 =  strcat(conditionname,num2str(degree),'deg_n1_RGB.mat');
                %file_Ref2 =  strcat(conditionname,num2str(degree),'deg_n1_RGB.mat');
                %file_Ref3 =  strcat(conditionname,num2str(degree),'deg_n1_RGB.mat');
                PM2(h).xCurrent
                file_Test =  strcat(pwd,'/Stimuli/Img_ang',num2str(degree),'_n',num2str(PM2(h).xCurrent),'.mat');
            
                %file_Ref =  strcat('UniformPatch_',num2str(degree),'deg_n1.mat');
                %file_Test =  strcat('UniformPatch_',num2str(degree),'deg_n',num2str(Intensity(h,Trial_cnt(h,1))),'.mat');
                
                factor = 3.2117*1.5/4.7;
                % 225.5 pixel 3 degree   
                % 188.9 pixel 2.5 degree   
            
                load(file_Test);Image(:,:,:,1) = stim;
                load(file_Ref1);Image(:,:,:,2) = stim;
            %     load(file_Ref2);Image(:,:,:,3) = stim;
            %     load(file_Ref3);Image(:,:,:,4) = stim;
            
                %file_TestName{Trial,h} = file_Test;
            
                imsize = size(Image);
            
                % Image to display stimuli.
                StimuliImage = zeros(imsize(2)*2+gapsize,imsize(1)*2+gapsize,3);
            
                arrangement = randperm(2);
                ansIndex2(h,Trial_cnt2(h,1)) = find(arrangement==1);
            
                % n
                % 1 3
                % 2 4
                r1 = 1:imsize(1);c1 = 1:imsize(2);
                r2 = imsize(1)+2*gapsize:2*imsize(1)+2*gapsize-1;c2 = 1:imsize(2);
                r3 = 1:imsize(1);c3 = imsize(2)+2*gapsize:2*imsize(2)+2*gapsize-1;
                r4 = imsize(1)+2*gapsize:2*imsize(1)+2*gapsize-1;c4 = imsize(1)+2*gapsize:2*imsize(1)+2*gapsize-1;
                
                StimuliImage(r2,c2,:) = Image(:,:,:,arrangement(1));
                StimuliImage(r4,c4,:) = Image(:,:,:,arrangement(2));
            %     StimuliImage(r3,c3,:) = Image(:,:,:,arrangement(3));
            %     StimuliImage(r4,c4,:) = Image(:,:,:,arrangement(4));
            
                StimuliImage = max(StimuliImage,0);
                StimuliImage = min(StimuliImage,1);
                Imsize = size(StimuliImage);
            
                crsSetDrawPage(CRS.VIDEOPAGE,1,1);
                crsSetCommand(CRS.CYCLEPAGEDISABLE);
                crsDrawMatrix42bitColour(coordsTarget,double(StimuliImage(:,1:2:Imsize(2),:)));
                crsSetZoneDisplayPage(CRS.VIDEOPAGE,1);
                
                pause(StimulusDuration);
                
                makeChromTemporalNoiseMask_TM(2,noise);
                crsResponseBoxFlush;
            
                button = -1;
                while button == -1;
                    BR = crsGetResponse(CRS.respCEDRUS);
                    pause(0.01);
            %         if (BR.Top == 1)||(BR.Left == 1)||(BR.Bottom == 1)||(BR.Right == 1)
            %             button = 1;
            %         end
                    if (BR.Left == 1)||(BR.Right == 1)
                        button = 1;
                    end
                pause(0.01);
                end
            
            %     if (BR.Top == 1) % top button
            %         response(h,Trial_cnt(h,1)) = 1;
            %     end
                if (BR.Left == 1) % left button
                    response(h,Trial_cnt2(h,1)) = 2;
                end
            %     if (BR.Bottom == 1) % bottom button
            %         response(h,Trial_cnt(h,1)) = 4;
            %     end
                if (BR.Right == 1) % right button
                    response(h,Trial_cnt2(h,1)) = 3;
                end
            
                %calculate performance or detect forced exit
                if response(h,Trial_cnt2(h,1)) == ansIndex(h,Trial_cnt2(h,1))
                    sound(s15,44100);
                    if practice == 0
                        correct2(h,Trial_cnt2(h,1)) = 1;
                    end
                elseif response(h,Trial_cnt2(h,1)) ~= ansIndex(h,Trial_cnt2(h,1))
                    sound(sx,44100);
                    if practice == 0
                        correct2(h,Trial_cnt2(h,1)) = 0;
                    end
                end
                pause(PauseAfterTrial);
            
                % Update intensity based on palamedes
                if practice == 0
                    PM2(h) = PAL_AMPM_updatePM(PM2(h), correct2(h,Trial_cnt2(h,1))); %update PM structure
                    % end if last 5 responses are stable
                    if Trial_cnt2(h,1) > 10 && mean(PM2(h).seThreshold(Trial_cnt2(h,1)-5:Trial_cnt2(h,1)-1)) < 2
                        PM2(h).stop = 1;
                    end
                
                    %if (PM(1).stop)&&(PM(2).stop)&&(PM(3).stop)&&(PM(4).stop)&&(PM(5).stop)&&(PM(6).stop)&&(PM(7).stop)&&(PM(8).stop) == 1
                    
                    for i = 1:noh
                        stop2(i,1) = PM2(i).stop;
                    end
                    
                    if (sum(stop2) == noh && sum(stop) == noh)
                        pause(0.5);
                        [tada, Fs] = wavread('C:\Windows\Media\tada.wav');
                        wavplay(2*tada,Fs);
                        vsgBlack
                    end
                end
                
                Trial_cnt2(h,1) = Trial_cnt2(h,1)+1;
                
                if min(Trial_cnt2(:,1) == (nop+1)*ones(noh,1))==1 && practice == 1
                    Trial_cnt2(:,1) = ones(noh,1);
                    practice = 0;
                end
            end
        end
    end 
end
end %of all interleaved staircase
 %of experiment

ShowCursor;
disp('Finished experiment');

%% Top


for i = 1:noh
    Threshold(i,:) = MB_AllHue(i,nor-round(PM(i).threshold(Trial_cnt(h,1)-1))+1,:);   
    Threshold_Number(i,:) = PM(i).threshold(Trial_cnt(h,1)-1);   
    %Threshold(i,1) = PM(i).threshold(end);
    %Threshold(i,2) = PM(i).threshold(end);
    
    Slope(i,1) = 10.^PM(i).slope(Trial_cnt(h,1)-1);   %PM.slope is in log10 units of beta parameter
end

result.threshold = Threshold;
result.refnumber = Threshold_Number;
result.slope = Slope;

currentFolder = pwd;
cd('C:\Data\Experiments\Reflectance Discrimination\Results');
filename = strcat('Top_session',num2str(session),'_',observer,'.mat');
%PMOut = PM;

for i = 1:noh
    PMOut(i).response = PM(i).response;
    PMOut(i).pdf = PM(i).pdf;
    PMOut(i).threshold = PM(i).threshold;
    PMOut(i).slope = PM(i).slope;
    PMOut(i).guess = PM(i).guess;
    PMOut(i).lapse = PM(i).lapse;
    PMOut(i).seThreshold = PM(i).seThreshold;
    PMOut(i).seSlope = PM(i).seSlope;
    PMOut(i).seGuess = PM(i).seGuess;
    PMOut(i).seLapse = PM(i).seLapse;
    PMOut(i).x = PM(i).x;
    PMOut(i).numTrials = PM(i).numTrials;
end

save(filename,'result','PMOut')
cd(currentFolder);

f = 1;

for i = 1:noh
    degree = (i-1)*360/noh;
    if i == 1
        fig(f) = figure(f);f = f + 1;
    end
    subplot(2,4,i);
    t(i,1:length(PM(i).x)) = 1:length(PM(i).x);
    plot(t(i,1:length(PM(i).x)),PM(i).x,'k');
    hold on;
    plot(t(i,PM(i).response == 1),PM(i).x(PM(i).response == 1),'ko', 'MarkerFaceColor','k');
    plot(t(i,PM(i).response == 0),PM(i).x(PM(i).response == 0),'ko', 'MarkerFaceColor','r');
    
    plot(1:length(PM(i).x)-1,PM(i).threshold,'b');

    set(gca,'FontSize',16);
    axis([0 max(t(i,1:length(PM(i).x)))+1 min(PM(i).x)-(max(PM(i).x)-min(PM(i).x))/10 max(PM(i).x)+(max(PM(i).x)-min(PM(i).x))/10]);
    xlabel('Trial');
    ylabel('Stimulus Intensity');
end

for i = 1:noh
    if i == 1
        fig(f) = figure(f);f = f + 1;
    end
    subplot(2,4,i);
    image(PAL_Scale0to1(PM(i).pdf)*64);hold on;
    set(gca,'FontSize',16);
    %axis([0 grain 0 grain]);
    xlabel('beta');
    ylabel('alpha');
end

fig(f) = figure(f);f = f + 1;
scatter(Threshold(:,1),Threshold(:,2));hold on;
scatter(0.7078,1,'rx');hold on;



%% Bottom

for i = 1:noh
    Threshold2(i,:) = MB_AllHue(i,nor-round(PM2(i).threshold(Trial_cnt2(h,1)-1))+1,:);   
    Threshold_Number2(i,:) = PM2(i).threshold(Trial_cnt2(h,1)-1);   
    %Threshold(i,1) = PM(i).threshold(end);
    %Threshold(i,2) = PM(i).threshold(end);
    
    Slope2(i,1) = 10.^PM2(i).slope(Trial_cnt2(h,1)-1);   %PM.slope is in log10 units of beta parameter
end

result2.threshold = Threshold2;
result2.refnumber = Threshold_Number2;
result2.slope = Slope2;

currentFolder = pwd;
cd('C:\Data\Experiments\Reflectance Discrimination\Results');
filename = strcat('Bottom_session',num2str(session),'_',observer,'.mat');
%PMOut = PM;

for i = 1:noh
    PMOut2(i).response = PM2(i).response;
    PMOut2(i).pdf = PM2(i).pdf;
    PMOut2(i).threshold = PM2(i).threshold;
    PMOut2(i).slope = PM2(i).slope;
    PMOut2(i).guess = PM2(i).guess;
    PMOut2(i).lapse = PM2(i).lapse;
    PMOut2(i).seThreshold = PM2(i).seThreshold;
    PMOut2(i).seSlope = PM2(i).seSlope;
    PMOut2(i).seGuess = PM2(i).seGuess;
    PMOut2(i).seLapse = PM2(i).seLapse;
    PMOut2(i).x = PM2(i).x;
    PMOut2(i).numTrials = PM2(i).numTrials;
end

save(filename,'result','PMOut2')
cd(currentFolder);

f = 1;

for i = 1:noh
    degree = (i-1)*360/noh;
    if i == 1
        fig2(f) = figure(f);f = f + 1;
    end
    subplot(2,4,i);
    t(i,1:length(PM2(i).x)) = 1:length(PM2(i).x);
    plot(t(i,1:length(PM2(i).x)),PM2(i).x,'k');
    hold on;
    plot(t(i,PM2(i).response == 1),PM2(i).x(PM(i).response == 1),'ko', 'MarkerFaceColor','k');
    plot(t(i,PM2(i).response == 0),PM2(i).x(PM(i).response == 0),'ko', 'MarkerFaceColor','r');
    
    plot(1:length(PM2(i).x)-1,PM2(i).threshold,'b');

    set(gca,'FontSize',16);
    axis([0 max(t(i,1:length(PM2(i).x)))+1 min(PM2(i).x)-(max(PM2(i).x)-min(PM2(i).x))/10 max(PM2(i).x)+(max(PM2(i).x)-min(PM2(i).x))/10]);
    xlabel('Trial');
    ylabel('Stimulus Intensity');
end

for i = 1:noh
    if i == 1
        fig2(f) = figure(f);f = f + 1;
    end
    subplot(2,4,i);
    image(PAL_Scale0to1(PM2(i).pdf)*64);hold on;
    set(gca,'FontSize',16);
    %axis([0 grain 0 grain]);
    xlabel('beta');
    ylabel('alpha');
end

fig2(f) = figure(f);f = f + 1;
scatter(Threshold2(:,1),Threshold2(:,2));hold on;
scatter(0.7078,1,'rx');hold on;
