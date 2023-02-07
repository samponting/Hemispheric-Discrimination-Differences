%function data1 = ModelingRD_PlotDatawithResults(En,Glossiness,xlim,ylim,color)

% load various files
load([pwd,'/Associated Files/blackbodylocus100to1000000_MB.mat'])
load([pwd,'/Associated Files/Basisvectors.mat'])
%load([pwd,'/Associated Files/Results_HumanThreshold'])
%load([pwd,'/Associated Files/MBandRGB_En1to4'])

numLightprobes = 1:10;
%mode = 'Upper';% 'Lower' or 'Upper'
% git
for mode = {'Upper','Lower'}
    for lps = 1%:numLightprobes
        xRGB.(['Env',num2str(lps)]) = load([pwd,'/HyperspectralLightprobe/HS_En',num2str(lps),'.mat']);
        sRGB.(['Env',num2str(lps)]) = xRGB.(['Env',num2str(lps)])./max(xRGB.(['Env',num2str(lps)]),[],'all');
        RGBupper.(['Env',num2str(lps)]) = sRGB.(['Env',num2str(lps)])(1:floor(size(sRGB.(['Env',num2str(lps)]),1)/2),:,:);
        RGBlower.(['Env',num2str(lps)]) = sRGB.(['Env',num2str(lps)])(floor(size(sRGB.(['Env',num2str(lps)]),1)/2)+1:size(sRGB.(['Env',num2str(lps)]),1),:,:);
        RGBup = sRGB.(['Env',num2str(lps)])(1:floor(size(sRGB.(['Env',num2str(lps)]),1)/2),:,:);
        RGBlo = sRGB.(['Env',num2str(lps)])(floor(size(sRGB.(['Env',num2str(lps)]),1)/2)+1:size(sRGB.(['Env',num2str(lps)]),1),:,:);
        
        imwrite(RGBup,[pwd,'/Figures/Env',num2str(lps),'_Img_Upper.png'])
        imwrite(RGBlo,[pwd,'/Figures/Env',num2str(lps),'_Img_Lower.png'])
        if strcmp(mode,'Upper')
            MBUpper.(['Env',num2str(lps)]) = RGBtoMB(RGBupper.(['Env',num2str(lps)]));
            RGBUpper = RGBupper;
        elseif strcmp(mode,'Lower')
            MBLower.(['Env',num2str(lps)]) = RGBtoMB(RGBlower.(['Env',num2str(lps)]));
            RGBLower = RGBlower;
        end
    end
end

%Threshold.HumanAverage = getfield(Average_HumanThreshold,['Environment',num2str(En)],Glossiness);

%% Figure 

for env = 1:numLightprobes

    v_En1_temp = reshape(MBUpper.(['Env' num2str(env)]),size(MBUpper.(['Env' num2str(env)]),1)*size(MBUpper.(['Env' num2str(env)]),2),size(MBUpper.(['Env' num2str(env)]),3));
    v_En1 = DeleteZeroLum(v_En1_temp,0.2);

    fig = figure;
    
    a = polyfit(v_En1(:,1)-0.7078,v_En1(:,2)-1,1);
    v = DeleteZeroLum(v_En1_temp,0.05);
    EnMB = reshape(MBUpper.(['Env' num2str(env)]),size(MBUpper.(['Env' num2str(env)]),1)*size(MBUpper.(['Env' num2str(env)]),2),3);EnRGB = reshape(RGBUpper.(['Env' num2str(env)]),size(RGBUpper.(['Env' num2str(env)]),1)*size(RGBUpper.(['Env' num2str(env)]),2),3);
    En_temp = reshape(MBUpper.(['Env' num2str(env)]),size(MBUpper.(['Env' num2str(env)]),1)*size(MBUpper.(['Env' num2str(env)]),2),size(MBUpper.(['Env' num2str(env)]),3));
    v_En = DeleteZeroLum(En_temp,0.2);
    a = polyfit(v_En(:,1)-0.7078,v_En(:,2)-1,1);
    v = DeleteZeroLum(En_temp,0.05);
    En_mean = [mean(v(:,1)),mean(v(:,2))];
    %EnMB = reshape(MB.(['Env' num2str(env)]),size(MB.(['Env' num2str(env)]),1)*size(MB.(['Env' num2str(env)]),2),3);EnRGB = reshape(RGB.(['Env' num2str(env)]),size(RGB.(['Env' num2str(env)]),1)*size(RGB.(['Env' num2str(env)]),2),3);
%     scatter(EnMB(:,1),EnMB(:,2),5,EnRGB./max(EnRGB,[],2)*0.8,'o','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.5);hold on
%     
%     a(2) = 1-0.7078*a(1);
%     line([0 100],[a(2) a(1)*100+a(2)],'LineWidth',1.5,'Color',[1 0.2 0.2],'LineStyle','--');hold on
%     plot(bbl_MB(:,2),bbl_MB(:,3),'Color',[0 0 1],'LineWidth',1.5,'LineStyle','--')
%     scatter(0.7078,1,200,[0.3 0.3 0.3],'+','LineWidth',1.5);hold on
%     scatter(0.6891,1.0714,200,[0.3 0.3 0.3],'+','LineWidth',1.5);hold on
% 
%     scatter(En_mean(1),En_mean(2),230,[0.3 0.3 0.3],'x','LineWidth',1.5);hold on
%     [~,data] = fitellipse(Threshold.HumanAverage(:,1),Threshold.HumanAverage(:,2));
%     p = plot(data(:,1),data(:,2),'-');hold on;


%%

    v_En1_temp_2 = reshape(MBLower.(['Env' num2str(env)]),size(MBLower.(['Env' num2str(env)]),1)*size(MBLower.(['Env' num2str(env)]),2),size(MBLower.(['Env' num2str(env)]),3));
    v_En1_2 = DeleteZeroLum(v_En1_temp_2,0.2);

    %fig = figure;
    
    a_2 = polyfit(v_En1_2(:,1)-0.7078,v_En1_2(:,2)-1,1);
    v_2 = DeleteZeroLum(v_En1_temp_2,0.05);
    EnMB_2 = reshape(MBLower.(['Env' num2str(env)]),size(MBLower.(['Env' num2str(env)]),1)*size(MBLower.(['Env' num2str(env)]),2),3);EnRGB_2 = reshape(RGBLower.(['Env' num2str(env)]),size(RGBLower.(['Env' num2str(env)]),1)*size(RGBLower.(['Env' num2str(env)]),2),3);
    En_temp_2 = reshape(MBLower.(['Env' num2str(env)]),size(MBLower.(['Env' num2str(env)]),1)*size(MBLower.(['Env' num2str(env)]),2),size(MBLower.(['Env' num2str(env)]),3));
    v_En_2 = DeleteZeroLum(En_temp_2,0.2);
    a_2 = polyfit(v_En_2(:,1)-0.7078,v_En_2(:,2)-1,1);
    v_2 = DeleteZeroLum(En_temp_2,0.05);
    En_mean_2 = [mean(v_2(:,1)),mean(v_2(:,2))];
    %EnMB = reshape(MB.(['Env' num2str(env)]),size(MB.(['Env' num2str(env)]),1)*size(MB.(['Env' num2str(env)]),2),3);EnRGB = reshape(RGB.(['Env' num2str(env)]),size(RGB.(['Env' num2str(env)]),1)*size(RGB.(['Env' num2str(env)]),2),3);
%     scatter(EnMB_2(:,1),EnMB_2(:,2),5,EnRGB_2./max(EnRGB_2,[],2)*0.8,'o','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.5);hold on
%     
%     a_2(2) = 1-0.7078*a_2(1);
    %line([0 100],[a_2(2) a_2(1)*100+a_2(2)],'LineWidth',1.5,'Color',[1 0.2 0.2],'LineStyle','--');hold on
    %plot(bbl_MB(:,2),bbl_MB(:,3),'Color',[0 0 1],'LineWidth',1.5,'LineStyle','--')
    %scatter(0.7078,1,200,[0.3 0.3 0.3],'+','LineWidth',1.5);hold on
    %scatter(En_mean_2(1),En_mean_2(2),230,[0.3 0.3 0.3],'x','LineWidth',1.5);hold on
    %[~,data] = fitellipse(Threshold.HumanAverage(:,1),Threshold.HumanAverage(:,2));
    %p = plot(data(:,1),data(:,2),'-');hold on;
    fig.PaperType       = 'a4';
    fig.PaperUnits      = 'centimeters';
    fig.PaperPosition   = [0,10,8.85,8.45];
    fig.Units           = 'centimeters';
    fig.Position        = [0,10,8.85,8.45];
    fig.Color           = 'w';
    fig.InvertHardcopy  = 'off';
    % 
    ax = gca;
    ax.XTick = [0 1];ax.XLim = [0 1];
    ax.YTick = [0 9];ax.YLim = [0 9];
%     ax.XTick = [min([elptUpper(:,1)' elptLower(:,1)']) max([elptUpper(:,1)' elptLower(:,1)'])];ax.XLim = [min([elptUpper(:,1)' elptLower(:,1)']) max([elptUpper(:,1)' elptLower(:,1)'])];
%     ax.YTick = [min([elptUpper(:,2)' elptLower(:,2)']) max([elptUpper(:,2)' elptLower(:,2)'])];ax.YLim = [min([elptUpper(:,2)' elptLower(:,2)']) max([elptUpper(:,2)' elptLower(:,2)'])];
%     % 
%     ax.XTickLabel = char(num2str(round(min([elptUpper(:,1)' elptLower(:,1)'])*1000)/1000),num2str(round(max([elptUpper(:,1)' elptLower(:,1)'])*1000)/1000));
%     ax.YTickLabel = char(num2str(round(min([elptUpper(:,2)' elptLower(:,2)'])*100)/100),num2str(round(max([elptUpper(:,2)' elptLower(:,2)'])*100)/100));
%     ax.XTickLabel = [];
%     ax.YTickLabel = [];
%     ax.XTick = [];
%     ax.YTick = [];
    xlabel('L/(L+M)');
    ylabel('S/(L+M)');
    % 
    ax.FontName = 'Arial';
    ax.FontSize = 12;
    % 
    ax.LineWidth = 0.4;
    ax.Units = 'centimeters';
    ticklengthcm(ax,0.3)
    axis square;
    ax.Position = [2.1 1.85 6 6];
    
    box off

    %saveas(gcf,[pwd,'/Figures/Env',num2str(env),'_Lower.png'])
    
    MBBoundary = zeros([length(EnMB) 3]);
    x = 1;
    for i = 1:length(EnMB)
        if 0 <= EnMB(i,1)  && 0 <= EnMB(i,2)
            MBBoundary(x,:) = EnMB(i,:);
            x = x+1;
        end
    end
    MBBoundary = MBBoundary(1:x-1, :);
    EnMB = MBBoundary;


    MBBoundary2 = zeros([length(EnMB_2) 3]);
    x = 1;
    for i = 1:length(EnMB_2)
        if 0 <= EnMB_2(i,1)  && 0 <= EnMB_2(i,2)
            MBBoundary2(x,:) = EnMB_2(i,:);
            x = x+1;
        end
    end
    MBBoundary2 = MBBoundary2(1:x-1, :);
    EnMB_2 = MBBoundary2;
    %elptUpper=ellipsedata(cov(EnMB(:,1),EnMB(:,2)),[mean(EnMB(:,1)),mean(EnMB(:,2))],length(EnMB),1);
    elptUpper=ellipsedata(cov(EnMB(:,1),EnMB(:,2)),[mean(EnMB(:,1)),mean(EnMB(:,2))],8,1);

    p = plot(elptUpper(:,1),elptUpper(:,2),'-');hold on
    p(1).LineWidth = 2;p(1).Color = [1 0 1];
%     

    %elptLower=ellipsedata(cov(EnMB_2(:,1),EnMB_2(:,2)),[mean(EnMB_2(:,1)),mean(EnMB_2(:,2))],length(EnMB_2),1);
    elptLower=ellipsedata(cov(EnMB_2(:,1),EnMB_2(:,2)),[mean(EnMB_2(:,1)),mean(EnMB_2(:,2))],8,1);

    p = plot(elptLower(:,1),elptLower(:,2),'-');hold on

    p(1).LineWidth = 2;p(1).Color = [0 0 1];
    
    %saveas(gcf,[pwd,'/Figures/Env',num2str(env),'_BiVar_Lower.png'])

    save([pwd,'/Results/Bivariate_',num2str(env),'_Lower.mat'],"elptLower")
    save([pwd,'/Results/Bivariate_',num2str(env),'_Upper.mat'],"elptUpper")



    fig.PaperType       = 'a4';
    fig.PaperUnits      = 'centimeters';
    fig.PaperPosition   = [0,10,8.85,8.45];
    fig.Units           = 'centimeters';
    fig.Position        = [0,10,8.85,8.45];
    fig.Color           = 'w';
    fig.InvertHardcopy  = 'off';
    % 
    ax = gca;
    ax.XTick = [0 1];ax.XLim = [0 1];
    ax.YTick = [0 9];ax.YLim = [0 9];
%     ax.XTick = [min([elptUpper(:,1)' elptLower(:,1)']) max([elptUpper(:,1)' elptLower(:,1)'])];ax.XLim = [min([elptUpper(:,1)' elptLower(:,1)']) max([elptUpper(:,1)' elptLower(:,1)'])];
%     ax.YTick = [min([elptUpper(:,2)' elptLower(:,2)']) max([elptUpper(:,2)' elptLower(:,2)'])];ax.YLim = [min([elptUpper(:,2)' elptLower(:,2)']) max([elptUpper(:,2)' elptLower(:,2)'])];
%     % 
%     ax.XTickLabel = char(num2str(round(min([elptUpper(:,1)' elptLower(:,1)'])*1000)/1000),num2str(round(max([elptUpper(:,1)' elptLower(:,1)'])*1000)/1000));
%     ax.YTickLabel = char(num2str(round(min([elptUpper(:,2)' elptLower(:,2)'])*100)/100),num2str(round(max([elptUpper(:,2)' elptLower(:,2)'])*100)/100));
%     ax.XTickLabel = [];
%     ax.YTickLabel = [];
%     ax.XTick = [];
%     ax.YTick = [];
    xlabel('L/(L+M)');
    ylabel('S/(L+M)');
    % 
    ax.FontName = 'Arial';
    ax.FontSize = 12;
    % 
    ax.LineWidth = 0.4;
    ax.Units = 'centimeters';
    ticklengthcm(ax,0.3)
    axis square;
    ax.Position = [2.1 1.85 6 6];
    
    box off
%     saveas(gcf,[pwd,'/Figures/Env',num2str(env),'_Combined.png'])   

% comment 

end

function [MBout,Id] = DeleteZeroLum(MBin,l)
    Id = find(MBin(:,3)>l);
    MBout = MBin(Id,:);
end