global CRS
crsInit;
crsSetVideoMode(CRS.HYPERCOLOURMODE);
PR670init('COM12');
wls = 380:5:780;
on = ones([crsGetScreenHeightPixels crsGetScreenWidthPixels]);
off = zeros([crsGetScreenHeightPixels crsGetScreenWidthPixels]);
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
