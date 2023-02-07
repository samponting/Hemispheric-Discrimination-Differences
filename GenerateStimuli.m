clear
close all

angleListDeg = 0:45:360;
angleListDeg = angleListDeg(1:length(angleListDeg)-1);
angleListRad = deg2rad(angleListDeg);
startDist = 0.1;
%white = ones([1 1 3]);
%white_MB = squeeze(RGBtoMB(white));
white_MB = [0.6891, 1.0714, 1];
MBwp(1,1,1) = white_MB(1);
MBwp(1,1,2) = white_MB(2);
MBwp(1,1,3) = white_MB(3);
white_RGB = squeeze(MBtoRGB(MBwp));
startPoints = zeros([2, length(angleListRad)]);
steps = 30;
rg_threshold = 1%0.0071;
yb_threshold = 1%0.1325;
for x = 1:steps
    for z = 1:length(angleListRad)
        [startPoints(1,z),startPoints(2,z)] = pol2cart(angleListRad(z),(x-1)*startDist/steps);
    %     scatter(startPoints(1,z),startPoints(2,z));hold on
    %     ax = gca;
    %     %xlim([0 1])
    %     %ylim([0 1])
    %     plot([0,startPoints(1,z)],[0,startPoints(2,z)])
    end

    MBstartPoints(1,:) = startPoints(1,:)+white_MB(1)/rg_threshold;
    MBstartPoints(2,:) = startPoints(2,:)+white_MB(2)/yb_threshold;

    mat = ones(299);
    grayMat = ones(299);
    mat = gauss2d(mat,500,[149.5, 149.5]);
    % figure()
    % imshow(mat)
    % figure()
    RGB = [];
    for hueNum = 1:length(angleListRad)
        %fig = figure();
        hue(1,1,1) = MBstartPoints(1,hueNum)*rg_threshold;
        hue(1,1,2) = MBstartPoints(2,hueNum)*yb_threshold;
        hue(1,1,3) = 1;

        RGB(1:3,hueNum) = squeeze(MBtoRGB(hue));
        scatter(MBstartPoints(2,hueNum),MBstartPoints(1,hueNum),10,RGB(1:3,hueNum)','filled');hold on
        ax = gca;
        xlabel('L/(L+M)')
        ylabel('S/(L+M)')
        stim(:,:,1) = mat.*0.5*RGB(1,hueNum) + (1-mat).*0.5*white_RGB(1);
        stim(:,:,2) = mat.*0.5*RGB(2,hueNum) + (1-mat).*0.5*white_RGB(2);
        stim(:,:,3) = mat.*0.5*RGB(3,hueNum) + (1-mat).*0.5*white_RGB(3);
        %stim = stim/max(stim,[],'all');
        %stim = stim*100;
%         stim2(:,:,1) = (stim(:,:,1).*mat) + (grayMat.*(1-mat));
%         stim2(:,:,2) = (stim(:,:,2).*mat) + (grayMat.*(1-mat));
%         stim2(:,:,3) = (stim(:,:,3).*mat) + (grayMat.*(1-mat));
%         stim2 = stim2./2;
        
        %figure()
        %imshow(stim)
        %save([pwd,'/Stimuli/Img_ang',num2str(angleListDeg(hueNum)),'_n',num2str(x),'.mat'],"stim")
        %close(fig)
    end
end

% figure()
% for z = 1:length(angleListRad)
%     
%     scatter(MBstartPoints(2,z),MBstartPoints(1,z),'k');hold on
%     ax = gca;
%     %xlim([0 1])
%     %ylim([0 1])
%     plot([white_MB(2)/yb_threshold,MBstartPoints(2,z)],[white_MB(1)/rg_threshold,MBstartPoints(1,z)],'Color',RGB(:,z),'LineWidth',2);
%     %p.LineStyle("Color", RGB)
%     xlabel('S/(L+M)')
%     ylabel('L/(L+M)')
% end


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