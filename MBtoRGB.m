function RGBImage = MBtoRGB(MBImage)

% Convert MacLeod-Boynton image to RGB image
% Input: MacLeod boynton image (m*n*3)
%        r = L/(L+M), b = S(L+M), and Luminance
% Output: RGB image

Matrix_LMS2RGB = [0.101802832302184,-0.0989750571347019,0.0113987555572253;...
  -0.0125252940574878,0.0367498534242654,-0.0102593696770454;...
  -0.000747942202009832,-0.000719986694654401,0.0233634225182479];

Lw = 0.692839;
Mw = 0.349676;
Sw = 2.146879448901693;

imsize = size(MBImage);

Luminance = MBImage(:,:,3);
r = MBImage(:,:,1);
b = MBImage(:,:,2);

%Lw = 1;Mw=1;Sw=1;
LMS(:,:,1) = (Luminance.*r)/Lw;
LMS(:,:,3) = (Luminance.*b)/Sw;
LMS(:,:,2) = (Luminance-LMS(:,:,1)*Lw)/Mw;

RGB_vector = (Matrix_LMS2RGB*(reshape(LMS,imsize(1)*imsize(2),3))')';

RGBImage_temp = reshape(RGB_vector,imsize(1),imsize(2),3);
%RGBImage = RGBImage_temp;
RGBImage = max(RGBImage_temp/max(RGBImage_temp(:)),0).^(1/2.2);
end
