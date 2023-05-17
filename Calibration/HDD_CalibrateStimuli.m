videoInit;
stimulusLevel = 100;
stimSize = 199;
CF = csvread('linss10e_1.csv');
CF = CF(1:391,:)';
CF(isnan(CF)) = 0;
CF = CF(2:4,1:5:391)';
Lw = 0.692839;
Mw = 0.349676;
Sw = 2.146879448901693;
rg_threshold = 0.0071;
yb_threshold = 0.1325;
startDist = 0.03*100;

%white_MB = FindWhitepoint;white_MB
white_MB = [0.7078/rg_threshold 1/yb_threshold 40]; % equal energy white
white_RGB = MB2RGB(white_MB);

for degree = 0:45:315
    angRad = deg2rad(degree);
    [unityVector(1), unityVector(2)] = pol2cart(angRad,(stimulusLevel*startDist)/100);
    stimCoord = white_MB(1:2) + unityVector;
    stimCoord(3) = white_MB(3);
    stimCoord(1) = stimCoord(1)*rg_threshold;
    stimCoord(2) = stimCoord(2)*yb_threshold;
 
    stim_RGB = MB2RGB(stimCoord);
    r = ones([crsGetScreenHeightPixels crsGetScreenWidthPixels]).* stim_RGB(1);
    g = ones([crsGetScreenHeightPixels crsGetScreenWidthPixels]).* stim_RGB(2);
    b = ones([crsGetScreenHeightPixels crsGetScreenWidthPixels]).* stim_RGB(3);
    stimulus = cat(3,r,g,b);
    crsSetDrawPage(1);
    crsDrawMatrix42bitColour([0 0],stimulus);
    crsSetDisplayPage(1);
%     recording = PR670measspd;
%     recording = recording(3:end);
% 
%     Lrec = CF(:,1)' * recording;
%     Mrec = CF(:,2)' * recording;
%     Srec = CF(:,3)' * recording;
% 
%     MBrec = [0 0 0];
%     MBrec(3) = (Lw*Lrec)+(Mw*Mrec);
%     MBrec(1) = (Lw*Lrec)/MBrec(3);
%     MBrec(2) = (Sw*Srec)/MBrec(3);
%     MBrec(3) = MBrec(3).*683;
%     disp(num2str(degree));
%     disp(['Estimated LMS: ', num2str(stimCoord)])
%     disp(['Recorded LMS: ',num2str(MBrec)])
    pause(2)
end
