load([pwd,'/Associated Files/blackbodylocus100to1000000_MB.mat'])
load([pwd,'/Associated Files/Basisvectors.mat'])
CFs = csvread('linss10e_1.csv');
Lw = 0.692839;
Mw = 0.349676;
Sw = 2.146879448901693;%0.1085;
CF = CFs(11:10:311,:,:);
numlps = 10;
if ~exist([pwd,'/Figures/'],'dir')
    mkdir([pwd,'/Figures/'])
end

mode = 'Upper';
% 'Upper' - Plots the distribution for the upper hemisphere
% 'Lower' - Plots the distribution for the lower hemisphere
% 'Whole' - Plots the distribution for the whole image

biVar = false;

% true - Plots bivariate distribution to distribution

rg_threshold = 0.0071; % scaling factors
yb_threshold = 0.1325;

rescale = true;

if rescale
    bbl_MB(:,2) = bbl_MB(:,2)/rg_threshold;
    bbl_MB(:,3) = bbl_MB(:,3)/yb_threshold;
end

savePlots = true;

%% Get Whitepoint

% equal energy white transformed into MB space using the CVRL 10 degree
% weightings L,M and S.
EEW = ones([1,441]);
whiteLMS = EEW*CFs(:,2:4);

MaB(1) = whiteLMS(1)*Lw/(Lw*whiteLMS(1)+Mw*whiteLMS(2));
MaB(2) = whiteLMS(3)*Sw/(Lw*whiteLMS(1)+Mw*whiteLMS(2));

if rescale
    MaB(1) = MaB(1)/rg_threshold;
    MaB(2) = MaB(2)/yb_threshold;
end
%% HS

for lps = 10
    load([pwd,'/HyperspectralLightprobe/HS_En',num2str(lps),'.mat']);
    xRGB = double(imread([pwd,'/HyperspectralLightprobe/HS_En',num2str(lps),'.png']));
    sRGB = xRGB./max(xRGB,[],'all');

    if strcmp(mode,'Upper')
        image = image(1:floor(size(image,1)/2),:,:);
        sRGB = sRGB(1:floor(size(sRGB,1)/2),:,:);
    elseif strcmp(mode,'Lower')
        image = image(floor(size(image,1)/2)+1:size(image,1),:,:);
        sRGB = sRGB(floor(size(sRGB,1)/2)+1:size(sRGB,1),:,:);
    end

    RGB = reshape(sRGB,size(sRGB,1)*size(sRGB,2),3);
    L = zeros([size(image,1) size(image,2)]);
    M = zeros([size(image,1) size(image,2)]);
    S = zeros([size(image,1) size(image,2)]);

    for x = 1:size(image,2)
        for y = 1:size(image,1)
            L(y,x) = CF(:,2)'*squeeze(image(y,x,:));
            M(y,x) = CF(:,3)'*squeeze(image(y,x,:));
            S(y,x) = CF(:,4)'*squeeze(image(y,x,:));
        end
    end
    
    L = reshape(L,size(L,1)*size(L,2),1);
    M = reshape(M,size(M,1)*size(M,2),1);
    S = reshape(S,size(S,1)*size(S,2),1);
    
    MB = zeros([length(L) 3]);
    MB(:,3) = ((Lw*L) + (Mw*M));
    MB(:,1) = (Lw*L)./MB(:,3);
    MB(:,2) = (Sw*S)./MB(:,3);
    
    if rescale
        MB(:,1) = MB(:,1)/rg_threshold;
        MB(:,2) = MB(:,2)/yb_threshold;
    end

    v_En = DeleteZeroLum(MB,0.2);
    a = polyfit(v_En(:,1)-MaB(1),v_En(:,2)-MaB(2),1);
    a(2) = MaB(2)-(MaB(1)*a(1));
    v = DeleteZeroLum(MB,0.05);
    En_mean = [mean(v(:,1)),mean(v(:,2))];
    fig = figure();

    scatter(MB(:,1),MB(:,2),2,RGB./max(RGB,[],2)*0.8,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.5);hold on
    
    %line([0 103],[a(2) a(1)*103+a(2)],'LineWidth',1.5,'Color',[1 0.2 0.2],'LineStyle','--');hold on
    %plot(bbl_MB(:,2),bbl_MB(:,3),'Color',[0 0 1],'LineWidth',1.5,'LineStyle','--')
    scatter(MaB(1),MaB(2),200,[0.3 0.3 0.3],'+','LineWidth',1.5);hold on
    %scatter(En_mean(1),En_mean(2),230,[0.3 0.3 0.3],'x','LineWidth',1.5);hold on
    tag = 'noEllipse';

    if biVar
        elpt=ellipsedata(cov(v(:,1),v(:,2)),[mean(v(:,1)),mean(v(:,2))],length(v),1);
        p = plot(elpt(:,1),elpt(:,2),'-');hold on
        p(1).LineWidth = 2;p(1).Color = [0 0 0];
        save([pwd,'/Results/Final/Bivariate_',num2str(lps),'_',mode,'.mat'],"elpt")
        tag = '+BiVar';
    end

    fig.PaperType       = 'a4';
    fig.PaperUnits      = 'centimeters';
    fig.PaperPosition   = [0,10,8.85,8.45];
    fig.Units           = 'centimeters';
    fig.Position        = [0,10,8.85,8.45];
    fig.Color           = 'w';
    fig.InvertHardcopy  = 'off';
    ax = gca;
    xbound = [0.65 0.74];
    ybound = [0.3 2.5];

    if rescale
       xbound = [90 105];
       ybound = [2,20];
    end
  
    ax.XTick = xbound;ax.XLim = xbound;ax.XTickLabel = xbound;
    ax.YTick = ybound;ax.YLim = ybound;ax.YTickLabel = ybound;
    xlabel('L/(L+M) rescaled');
    ylabel('S/(L+M) rescaled');
    ax.FontName = 'Arial';
    ax.FontSize = 12;
    ax.LineWidth = 0.4;
    ax.Units = 'centimeters';
    ticklengthcm(ax,0.3)
    axis square;
    ax.Position = [2.1 1.85 6 6];
    box off

    if savePlots == true
        if rescale
          saveas(gcf,[pwd,'/Figures/Env',num2str(lps),'_',mode,'_Dist','_Rescaled',tag,'.png'])
        else
          saveas(gcf,[pwd,'/Figures/Env',num2str(lps),'_',mode,'_Dist',tag,'.png'])
        end
    end
end

function [MBout,Id] = DeleteZeroLum(MBin,l)
    Id = find(MBin(:,3)>l);
    MBout = MBin(Id,:);
end
