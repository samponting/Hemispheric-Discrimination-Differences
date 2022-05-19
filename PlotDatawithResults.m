%function data1 = ModelingRD_PlotDatawithResults(En,Glossiness,xlim,ylim,color)

% load various files
load([pwd,'/Associated Files/blackbodylocus100to1000000_MB.mat'])
load([pwd,'/Associated Files/Basisvectors.mat'])
%load([pwd,'/Associated Files/Results_HumanThreshold'])
%load([pwd,'/Associated Files/MBandRGB_En1to4'])

numLightprobes = 10;
mode = 'Upper';% 'Lower' or 'Upper'


for lps = 1:numLightprobes
    xRGB.(['Env',num2str(lps)]) = double(imread([pwd,'/HyperspectralLightprobe/HS_En',num2str(lps),'.png']));
    sRGB.(['Env',num2str(lps)]) = xRGB.(['Env',num2str(lps)])./max(xRGB.(['Env',num2str(lps)]),[],'all');
    RGBupper.(['Env',num2str(lps)]) = sRGB.(['Env',num2str(lps)])(1:floor(size(sRGB.(['Env',num2str(lps)]),1)/2),:,:);
    RGBlower.(['Env',num2str(lps)]) = sRGB.(['Env',num2str(lps)])(floor(size(sRGB.(['Env',num2str(lps)]),1)/2)+1:size(sRGB.(['Env',num2str(lps)]),1),:,:);
    
    if strcmp(mode,'Upper')
        MB.(['Env',num2str(lps)]) = RGBtoMB(RGBupper.(['Env',num2str(lps)]));
        RGB = RGBupper;
    elseif strcmp(mode,'Lower')
        MB.(['Env',num2str(lps)]) = RGBtoMB(RGBlower.(['Env',num2str(lps)]));
        RGB = RGBlower;
    end
end


%Threshold.HumanAverage = getfield(Average_HumanThreshold,['Environment',num2str(En)],Glossiness);

%% Figure 

for env = 1:numLightprobes
    
    v_En1_temp = reshape(MB.(['Env' num2str(env)]),size(MB.(['Env' num2str(env)]),1)*size(MB.(['Env' num2str(env)]),2),size(MB.(['Env' num2str(env)]),3));
    %v_En2_temp = reshape(MB.Env2,size(MB.Env2,1)*size(MB_En2,2),size(MB_En2,3));
    % v_En3_temp = reshape(MB_En1,size(MB_En1,1)*size(MB_En1,2),size(MB_En1,3));
    % v_En3_temp(:,1) = 2*0.7078 - v_En1_temp(:,1);
    % v_En3_temp(:,2) = v_En1_temp(:,2);
    
    v_En1 = DeleteZeroLum(v_En1_temp,0.2);
    % v_En2 = DeleteZeroLum(v_En2_temp,0.2);
    
    % v_En3(:,1) = 2*0.7078 - v_En1(:,1);
    % v_En3(:,2) = v_En1(:,2);
    
    fig = figure;
    
    a = polyfit(v_En1(:,1)-0.7078,v_En1(:,2)-1,1);
    v = DeleteZeroLum(v_En1_temp,0.05);
    
    %En_mean = [mean(v(:,1)),mean(v(:,2))];
    
    
    EnMB = reshape(MB.(['Env' num2str(env)]),size(MB.(['Env' num2str(env)]),1)*size(MB.(['Env' num2str(env)]),2),3);EnRGB = reshape(RGB.(['Env' num2str(env)]),size(RGB.(['Env' num2str(env)]),1)*size(RGB.(['Env' num2str(env)]),2),3);
    
    En_temp = reshape(MB.(['Env' num2str(env)]),size(MB.(['Env' num2str(env)]),1)*size(MB.(['Env' num2str(env)]),2),size(MB.(['Env' num2str(env)]),3));
    v_En = DeleteZeroLum(En_temp,0.2);
    a = polyfit(v_En(:,1)-0.7078,v_En(:,2)-1,1);
    v = DeleteZeroLum(En_temp,0.05);
    En_mean = [mean(v(:,1)),mean(v(:,2))];
    EnMB = reshape(MB.(['Env' num2str(env)]),size(MB.(['Env' num2str(env)]),1)*size(MB.(['Env' num2str(env)]),2),3);EnRGB = reshape(RGB.(['Env' num2str(env)]),size(RGB.(['Env' num2str(env)]),1)*size(RGB.(['Env' num2str(env)]),2),3);
    % EnMB = MB.(['Env',num2str(En)]);
    % EnRGB = RGB.(['Env',num2str(En)]);
    
    scatter(EnMB(:,1),EnMB(:,2),5,EnRGB./max(EnRGB,[],2)*0.8,'o','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.5);hold on
    
    a(2) = 1-0.7078*a(1);
    line([0 100],[a(2) a(1)*100+a(2)],'LineWidth',1.5,'Color',[1 0.2 0.2],'LineStyle','--');hold on
    plot(bbl_MB(:,2),bbl_MB(:,3),'Color',[0 0 1],'LineWidth',1.5,'LineStyle','--')
    scatter(0.7078,1,200,[0.3 0.3 0.3],'+','LineWidth',1.5);hold on
    scatter(En_mean(1),En_mean(2),230,[0.3 0.3 0.3],'x','LineWidth',1.5);hold on
    %[~,data] = fitellipse(Threshold.HumanAverage(:,1),Threshold.HumanAverage(:,2));
    %p = plot(data(:,1),data(:,2),'-');hold on;

    elpt=ellipsedata(cov(EnMB(:,1),EnMB(:,2)),[mean(EnMB(:,1)),mean(EnMB(:,2))],length(EnMB),0.5);
    p = plot(elpt(:,1),elpt(:,2),'-');hold on
    p(1).LineWidth = 2;p(1).Color = [0.3 0.3 0.3];
    

    fig.PaperType       = 'a4';
    fig.PaperUnits      = 'centimeters';
    fig.PaperPosition   = [0,10,8.85,8.45];
    fig.Units           = 'centimeters';
    fig.Position        = [0,10,8.85,8.45];
    fig.Color           = 'w';
    fig.InvertHardcopy  = 'off';
    % 
    ax = gca;
    % ax.XTick = [0.695,0.7075,0.720];ax.XLim = [0 1];
    ax.XTick = [min(elpt(:,1)) max(elpt(:,1))];ax.XLim = [min(elpt(:,1)) max(elpt(:,1))];
    ax.YTick = [min(elpt(:,2)) max(elpt(:,2))];ax.YLim = [min(elpt(:,2)) max(elpt(:,2))];
    % 
    ax.XTickLabel = char(num2str(round(min([0.4 0.9])*1000)/1000),num2str(round(max([0.4 0.9])*1000)/1000));
    ax.YTickLabel = char(num2str(round(min([0.7 2.5])*100)/100),num2str(round(max([0.7 2.5])*100)/100));
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
    saveas(gcf,[pwd,'/Figures/Env',num2str(env),mode,'.png'])   



end

function [MBout,Id] = DeleteZeroLum(MBin,l)
    Id = find(MBin(:,3)>l);
    MBout = MBin(Id,:);
end