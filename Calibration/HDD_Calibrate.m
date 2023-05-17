%HDD Calibration

global CRS
crsInit;
crsSetVideoMode(CRS.HYPERCOLOURMODE);
PR670init('COM12');
wls = 380:5:780;
on = ones([crsGetScreenHeightPixels crsGetScreenWidthPixels]);
off = zeros([crsGetScreenHeightPixels crsGetScreenWidthPixels]);
steps = 0:0.03:1;
CF = csvread('linss10e_1.csv');
load('linCIE2008v10e_1.csv');
d65 = load('D65.txt');
Lw = 0.692839;
Mw = 0.349676;
Sw = 2.146879448901693;

%% Get D65

CF = CF(1:391,:)';
CF(isnan(CF)) = 0;
CF = CF(2:4,:);
d65 = d65(91:481,:);
d65(:,2) = d65(:,2)/max(d65(:,2));
d65 = d65(:,2)';
d65 = d65./2;

whiteSpd = ones(1,391).*0.3;

LMS = d65 * CF';
LMS = whiteSpd * CF';

MB(3) = (LMS(1)*Lw) + (LMS(2)*Mw);
MB(1) = (LMS(1)*Lw)/MB(3);
MB(2) = (LMS(3)*Sw)/MB(3);

%% Get Display++Primaries

crsSetDrawPage(CRS.VIDEOPAGE,1);
crsDrawMatrix42bitColour([0 0], cat(3,on,off,off));
crsSetZoneDisplayPage(CRS.VIDEOPAGE,1);
for i = 1:5
    Red(:,i) = PR670measspd;
end
aveRed = mean(Red,2);
crsDrawMatrix42bitColour([0 0], cat(3,off,on,off));
crsSetZoneDisplayPage(CRS.VIDEOPAGE,1);
for i = 1:5
    Green = PR670measspd;
end
aveGreen = mean(Green,2);
crsDrawMatrix42bitColour([0 0], cat(3,off,off,on));
crsSetZoneDisplayPage(CRS.VIDEOPAGE,1);
for i = 1:5
    Blue = PR670measspd;
end
aveBlue = mean(Blue,2);
Primaries = [wls' aveRed aveGreen aveBlue];
save('Display++Primaries.mat','Primaries')

%% Get Gamma


for i = 1:length(steps)
    crsDrawMatrix42bitColour([0 0], cat(3,on*steps(i),off,off));
    crsSetZoneDisplayPage(CRS.VIDEOPAGE,1);
    Red(:,i) = PR670measspd;
    RedLum(i) = Red(3:end,i)'*linCIE2008v10e_1(1:5:391,2);
    disp(['measuring red level',num2str(i),' /100'])
    scatter(i,Red(3:end,i)'*linCIE2008v10e_1(1:5:391,2),'filled');hold on
end
RedFit = fit(steps',RedLum','b*x^c+d');
Rcoef = coeffvalues(RedFit);Rcoef = Rcoef(2);
close all
for i = 1:length(steps)
    crsDrawMatrix42bitColour([0 0], cat(3,off,on*steps(i),off));
    crsSetZoneDisplayPage(CRS.VIDEOPAGE,1);
    Green(:,i) = PR670measspd;
    GreenLum(i) = Green(3:end,i)'*linCIE2008v10e_1(1:5:391,2);
    disp(['measuring green level',num2str(i),' /100'])
    scatter(i,Green(3:end,i)'*linCIE2008v10e_1(1:5:391,2),'filled');hold on
end
GreenFit = fit(steps',GreenLum','b*x^c+d');
Gcoef = coeffvalues(GreenFit);Gcoef = Gcoef(2);
close all
for i = 1:length(steps)
    crsDrawMatrix42bitColour([0 0], cat(3,off,off,on*steps(i)));
    crsSetZoneDisplayPage(CRS.VIDEOPAGE,1);
    Blue(:,i) = PR670measspd;
    BlueLum(i) = Blue(3:end,i)'*linCIE2008v10e_1(1:5:391,2);
    disp(['measuring blue level',num2str(i),' /100'])
    scatter(i,Blue(3:end,i)'*linCIE2008v10e_1(1:5:391,2),'filled');hold on
end
BlueFit = fit(steps',BlueLum','b*x^c+d');
Bcoef = coeffvalues(BlueFit);Bcoef = Bcoef(2);

%% MB to RGB
load('Display++Primaries.mat');
CF = CF(:,1:5:391)';
Primaries = Primaries(3:end,:);

Lr = Primaries(:,2)' * CF(:,1);
Lg = Primaries(:,3)' * CF(:,1);
Lb = Primaries(:,4)' * CF(:,1);

Mr = Primaries(:,2)' * CF(:,2);
Mg = Primaries(:,3)' * CF(:,2);
Mb = Primaries(:,4)' * CF(:,2);

Sr = Primaries(:,2)' * CF(:,3);
Sg = Primaries(:,3)' * CF(:,3);
Sb = Primaries(:,4)' * CF(:,3);

RGB2LMS_Matrix = [Lr, Lg, Lb;
    Mr, Mg, Mb;
    Sr, Sg, Sb];

RGB2LMS_Matrix = RGB2LMS_Matrix.*683;

LMS2RGB_Matrix = inv(RGB2LMS_Matrix);

RGB = LMS2RGB_Matrix * LMS';

% RGB(1) = RGB(1).^(1/Rcoef);
% RGB(2) = RGB(2).^(1/Gcoef);
% RGB(3) = RGB(3).^(1/Bcoef);

%% Check White
crsSetDrawPage(1)
crsDrawMatrix42bitColour([0 0], cat(3,RGB(1)*on,RGB(2)*on,RGB(3)*on));
crsSetDisplayPage(1);

recording = PR670measspd;

recording = recording(3:end);

Lrec = CF(:,1)' * recording;
Mrec = CF(:,2)' * recording;
Srec = CF(:,3)' * recording;

MBrec = [0 0 0];
MBrec(3) = (Lw*Lrec)+(Mw*Mrec);
MBrec(1) = (Lw*Lrec)/MBrec(3);
MBrec(2) = (Sw*Srec)/MBrec(3);


disp('Estimated LMS: ')
LMS
disp('Recorded LMS: ')
[Lrec,Mrec,Srec].*683




