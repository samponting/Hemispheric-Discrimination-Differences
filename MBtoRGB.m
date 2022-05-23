function RGBImage = ModelingRD_MBtoRGBImage(MBImage)

% Convert MacLeod-Boynton image to RGB image
% Input: MacLeod boynton image (m*n*3)
%        r = L/(L+M), b = S(L+M), and Luminance
% Output: RGB image

Matrix_LMS2RGB = [0.0902678377277668,-0.0884125119776831,0.00658074156160040;...
    -0.0122358781709693,0.0358560339886004,-0.00693733149151792;...
    -0.000520348637699630,-0.00110419888293445,0.0212184883479314];

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
RGBImage = max(RGBImage_temp/max(RGBImage_temp(:)),0).^(1/2.2);
end
