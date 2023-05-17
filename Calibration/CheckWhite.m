% Check White
global CRS
crsInit;
crsSetVideoMode(CRS.HYPERCOLOURMODE);
CF = csvread('linss10e_1.csv');
CF = CF(1:5:391,:);
CF(isnan(CF)) = 0;
on = ones([crsGetScreenHeightPixels crsGetScreenWidthPixels]);
white = white_RGB;
crsSetDrawPage(CRS.VIDEOPAGE,1);
crsDrawMatrix42bitColour([0 0], cat(3,white(1)*on,white(2)*on,white(3)*on));
crsSetZoneDisplayPage(CRS.VIDEOPAGE,1);
Lw = 0.692839;
Mw = 0.349676;
Sw = 2.146879448901693;

PR670init('COM11');
recording = PR670measspd;

recording = recording(3:end);

L = CF(:,2)' * recording;
M = CF(:,3)' * recording;
S = CF(:,4)' * recording;

MB = [0 0 0];
MB(3) = (Lw*L)+(Mw*M);
MB(1) = (Lw*L)/MB(3);
MB(2) = (Sw*S)/MB(3);

