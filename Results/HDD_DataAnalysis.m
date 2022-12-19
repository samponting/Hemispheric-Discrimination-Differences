%% Initialise
% load files and other resources
files = dir('*.mat'); % results files - make sure all desired files are 
% in the same folder as this file. Any excluded data should be in a
% different folder.
CFs = csvread('linss10e_1.csv'); % cone fundamentals (for whitepoint). 
% these are found on cvrl.com, look for the linear energy cone fundamentals
% with a 10 degree field size. You should have this file already for the
% lightprobe analysis.
load([pwd, '/bbl/blackbodylocus100to1000000_MB.mat']); % blackbody locus
% Takuma gave me this file. I'll be sending it to you along with this code
% - make sure it is in a folder within the directory in order to
% separate it from the actual results

% search through the files extracting each of the participant's names. 
pptnames = {};
nmCount = 1;
for i = 1:length(files)
    if ~any(strcmp(pptnames,extractBefore(files(i).name,'_')))
        pptnames{nmCount} = extractBefore(files(i).name,'_');
        nmCount = nmCount+1;
    end
end

% logical flag for saving figures
saveFig = 0;

%% Get Whitepoint

% equal energy white transformed into MB space using the CVRL 10 degree
% weightings L,M and S.
EEW = ones([1,441]);
whiteLMS = EEW*CFs(:,2:4);
Lw = 0.692839;
Mw = 0.349676;
Sw = 2.146879448901693;

MB(1) = whiteLMS(1)*Lw/(Lw*whiteLMS(1)+Mw*whiteLMS(2));
MB(2) = whiteLMS(3)*Sw/(Lw*whiteLMS(1)+Mw*whiteLMS(2));

%% by ppt
% by participant figures and results
fig = figure();
fig.Position = [0 0 1500 2000];
% 'tiledlayout' and 'nexttile' create a large figure with multiple subplots
% for each of the participants. we loop through the participants ('i') and
% iterate onto the next subplot using 'nexttile'.
t = tiledlayout('flow','TileSpacing','compact');
for i = 1:length(pptnames)
    nexttile
    pptFiles = dir([pptnames{i},'_*.mat']);
    if length(pptFiles) > 1
        for y = 1:length(pptFiles)
            load(pptFiles(y).name)
            thresholds(:,:,y) = result.threshold(:,1:2);
        end
    end

    % calculate the mean thresholds for each participant, then plot the
    % upper (1:8) and lower (9:16) thresholds.
    meanThreshold = mean(thresholds,3);
    scatter(meanThreshold(1:8,1),meanThreshold(1:8,2),25,[1 0 0.3],'filled','MarkerEdgeColor','none');hold on
    scatter(meanThreshold(9:16,1),meanThreshold(9:16,2),25,[0.3 0 1],'filled','MarkerEdgeColor','none');
    % plot EEW
    scatter(MB(1),MB(2),200,[0 0 0],'x')
    % fit ellipses and plot for each set of thresholds.
    [a,b] = fitellipse(meanThreshold(1:8,1),meanThreshold(1:8,2));
    plot(b(:,1),b(:,2),'-','Color',[1 0 0.3],'LineWidth',1.3)
    [a,b] = fitellipse(meanThreshold(9:16,1),meanThreshold(9:16,2));
    plot(b(:,1),b(:,2),'-','Color',[0.3 0 1],'LineWidth',1.3)
    % plot blackbody locus
    plot(bbl_MB(:,2),bbl_MB(:,3),'--','Color',[0 0 0]);
    % subplot formatting
    xlim([0.5 0.89])
    ylim([0.84 1.23])
    xlabel('L/(L+M)')
    ylabel('S/(L+M)')
    axis square
    ax = gca;
    ax.LineWidth = 2;
    ax.FontName = 'Ariel';
    ax.FontSize = 14;
    title(sprintf('Discrimination Ellipses: %s',pptnames{i}))
end

%% overall
% overall data plot, using the same method as the previous section
for i = 1:length(files)
    load(files(i).name)
    thresholds(:,:,i) = result.threshold(:,1:2);
end
nexttile

meanThreshold = mean(thresholds,3);
scatter(meanThreshold(1:8,1),meanThreshold(1:8,2),25,[1 0 0.3],'filled','MarkerEdgeColor','none');hold on
scatter(meanThreshold(9:16,1),meanThreshold(9:16,2),25,[0.3 0 1],'filled','MarkerEdgeColor','none');hold on
scatter(MB(1),MB(2),200,[0 0 0],'x')
[a,b] = fitellipse(meanThreshold(1:8,1),meanThreshold(1:8,2));
plot(b(:,1),b(:,2),'-','Color',[1 0 0.3],'LineWidth',1.3)
[a,b] = fitellipse(meanThreshold(9:16,1),meanThreshold(9:16,2));
plot(b(:,1),b(:,2),'-','Color',[0.3 0 1],'LineWidth',1.3)
plot(bbl_MB(:,2),bbl_MB(:,3),'--','Color',[0 0 0]);
lgd = legend({'Upper','Lower','Equal Energy White','','','Black Body Locus'}, ...
    'FontSize',12);
legend('boxoff')
lgd.Layout.Tile = 'east';
xlim([0.5 0.89])
ylim([0.84 1.23])
xlabel('L/(L+M)')
ylabel('S/(L+M)')
axis square
ax = gca;
ax.LineWidth = 2;
ax.FontName = 'Ariel';
ax.FontSize = 14;
title('Discrimination Ellipses: Mean')
if saveFig == 1
    saveas(fig,'DiscriminationEllipses.png')
end

%% STATS: AREA

fig = figure();
fig.Position = [0 0 1500 2000];

subplot(1,3,1)
for i = 1:size(thresholds,3)
    UpperEllipse(:,i) = fitellipse(thresholds(1:8,1,i),thresholds(1:8,2,i));
    LowerEllipse(:,i) = fitellipse(thresholds(9:16,1,i),thresholds(9:16,2,i));
end

% Mann-Whitney U-Test (nonparametric)
areaStats = mwwtest(UpperEllipse(6,:),LowerEllipse(6,:));

% bar graph showing area
b = bar(1:2,[mean(UpperEllipse(6,:)), mean(LowerEllipse(6,:))],0.7,'FaceColor',[0.3 0 1]);
set(gca,'XTickLabel',["Upper","Lower"]);hold on
errorbar([mean(UpperEllipse(6,:)), mean(LowerEllipse(6,:))],[std(UpperEllipse(6,:))./sqrt(length(UpperEllipse)), std(LowerEllipse(6,:))./sqrt(length(LowerEllipse))],'LineStyle','none','Color',[0 0 0],'LineWidth',1.5)
axis square
ax = gca;
ax.LineWidth = 2;
ax.FontName = 'Ariel';
ax.FontSize = 14;
ylim([0.0025 0.0042])
title('Area of Ellipses')
ylabel('Area')
legend(sprintf('p = %s',string(round(areaStats.p(1),4))),'Standard Error')
if saveFig == 1
    saveas(fig,'EllipseAreas.png')
end

%% STATS: ECCENTRICITY
% this section essentially repeats the above section, but for the
% eccentricity of the ellipses

subplot(1,3,2)
for i = 1:length(UpperEllipse)
    % the equation for eccentricity is sqrt(1-b^2/a^2) where a and b are
    % the minor and major axes. 
    eccentricityUpper(i) = sqrt(1-(min(UpperEllipse(3:4,i)).^2./max(UpperEllipse(3:4,i)).^2));
    eccentricityLower(i) = sqrt(1-(min(LowerEllipse(3:4,i)).^2./max(LowerEllipse(3:4,i)).^2));
end




eccentricityStats = mwwtest(eccentricityUpper,eccentricityLower);

b = bar(1:2,[mean(eccentricityUpper), mean(eccentricityLower)],0.7,'FaceColor',[0.3 0 1]);
set(gca,'XTickLabel',["Upper","Lower"]);hold on
errorbar([mean(eccentricityUpper), mean(eccentricityLower)],[std(eccentricityUpper)/sqrt(length(UpperEllipse)), std(eccentricityLower)/sqrt(length(LowerEllipse))],'LineStyle','none','Color',[0 0 0],'LineWidth',1.5)
axis square
ax = gca;
ax.LineWidth = 2;
ax.FontName = 'Ariel';
ax.FontSize = 14;
title('Eccentricity of Ellipses')
ylabel('Eccentricity')
legend(sprintf('p = %s',string(round(eccentricityStats.p(1),4))),'Standard Error')
ylim([0.997 1])
if saveFig == 1
    saveas(fig,'EllipseAreas.png')
end

%% Find Angles

% this section finds the orientation of the ellipses, as I found that the
% in-built angle from fitellipse.m was inconsistent. the angle found is the
% angle of elevation of the major axis from the positive x axis.

for i = 1:length(files)
    load(files(i).name)
    thresholds = result.threshold(:,1:2);

    [a,b] = fitellipse(thresholds(1:8,1),thresholds(1:8,2));
    dists = pdist2(a(1:2),b);
    [~,maxDistInd] = max(dists);
    maxCoord = b(maxDistInd,:);
    if maxCoord(2) < a(2)
        maxCoord = -(maxCoord - a(1:2))+a(1:2);
    end
    slope = (maxCoord(2)-a(2)) / (maxCoord(1)-a(1));
    angleUpper(i) = atand(slope);
    if angleUpper(i) < 0
        angleUpper(i) = angleUpper(i) + 180;
    end

    [a,b] = fitellipse(thresholds(9:16,1),thresholds(9:16,2));
    dists = pdist2(a(1:2),b);
    [~,maxDistInd] = max(dists);
    maxCoord = b(maxDistInd,:);
    if maxCoord(2) < a(2)
        maxCoord = -(maxCoord - a(1:2))+a(1:2);
    end
    slope = (maxCoord(2)-a(2)) / (maxCoord(1)-a(1));
    angleLower(i) = atand(slope);
    if angleLower(i) < 0
        angleLower(i) = angleLower(i) + 180;
    end
end

%% STATS: ORIENTATION
% repeats above sections for the orientation data.

subplot(1,3,3)

orientationStats = mwwtest(angleUpper,angleLower);

b = bar(1:2,[mean(angleUpper), mean(angleLower)],0.7,'FaceColor',[0.3 0 1]);
set(gca,'XTickLabel',["Upper","Lower"]);hold on
errorbar([mean(angleUpper), mean(angleLower)],[std(angleUpper)./sqrt(length(UpperEllipse)), std(angleLower)./sqrt(length(LowerEllipse))],'LineStyle','none','Color',[0 0 0],'LineWidth',1.5)
axis square
ax = gca;
ax.LineWidth = 2;
ax.FontName = 'Ariel';
ax.FontSize = 14;
title('Orientation of Ellipses')
ylabel('Orientation (Degrees)')
legend(sprintf('p = %s',string(round(orientationStats.p(1),4))),'Standard Error')
ylim([88 92])
if saveFig == 1
    saveas(fig,'EllipseAreas.png')
end
