angleListDeg = 0:45:360;
angleListDeg = angleListDeg(1:length(angleListDeg)-1);
angleListRad = deg2rad(angleListDeg);
startDist = 0.10;
white = [0.5 0.5 0.5];
white_MB = RGBtoMB(white);
startPoints = zeros([2, length(angleListRad)]);
for z = 1:length(angleListRad)
    [startPoints(1,z),startPoints(2,z)] = pol2cart(angleListRad(z),startDist);
    
%     scatter(startPoints(1,z),startPoints(2,z));hold on
%     ax = gca;
%     %xlim([0 1])
%     %ylim([0 1])
%     plot([0,startPoints(1,z)],[0,startPoints(2,z)])
end

MBstartPoints(1,:) = startPoints(1,:) + white_MB(1);
MBstartPoints(2,:) = startPoints(2,:) + white_MB(2);

mat = ones(1000);

mat = gauss2d(mat,50000,[500, 500]);
% figure()
% imshow(mat)
figure()
RGB = [];
for hueNum = 1:length(angleListRad)
    hue(1,1,1) = MBstartPoints(1,hueNum);
    hue(1,1,2) = MBstartPoints(2,hueNum);
    hue(1,1,3) = 0.5;

    RGB(1:3,hueNum) = squeeze(MBtoRGB(hue));
    
    stim(:,:,1) = mat.*RGB(1);
    stim(:,:,2) = mat.*RGB(2);
    stim(:,:,3) = mat.*RGB(3);
    %figure()
    %imshow(stim)
    scatter(MBstartPoints(1,z),MBstartPoints(2,z));hold on
    ax = gca;
    %xlim([0 1])
    %ylim([0 1])
    plot([white_MB(1),MBstartPoints(1,hueNum)],[white_MB(2),MBstartPoints(2,hueNum)],'Color',RGB(:,hueNum));
    %p.LineStyle("Color", RGB)

end

for z = 1:length(angleListRad)
    
end


function mat = gauss2d(mat, sigma, center)
gsize = size(mat);
[R,C] = ndgrid(1:gsize(1), 1:gsize(2));
mat = gaussC(R,C, sigma, center);
end

function val = gaussC(x, y, sigma, center)
xc = center(1);
yc = center(2);
exponent = ((x-xc).^2 + (y-yc).^2)./(2*sigma);
val       = (exp(-exponent));
end