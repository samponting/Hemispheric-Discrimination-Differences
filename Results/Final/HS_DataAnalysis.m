%% Initialise
% Upperfiles = dir('*Upper.mat');
% Lowerfiles = dir('*Lower.mat');

% Logical flags
saveFig = 1;

%% STATS: AREA


fig = figure();
fig.Position = [0 0 150 200];

subplot(1,3,1)

for i = 1:10
    load(['Bivariate_' num2str(i) '_Upper.mat']);
    UpperEllipse(:,i) = fitellipse(elpt(:,1),elpt(:,2));
    load(['Bivariate_' num2str(i) '_Lower.mat']);
    LowerEllipse(:,i) = fitellipse(elpt(:,1),elpt(:,2));
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

ylim([0 5])
title('Area of Ellipses')
ylabel('Area')
legend(sprintf('p = %s',string(round(areaStats.p(1),4))),'Standard Error')
%if saveFig == 1
   %saveas(fig,'HSEllipseStats.png')
%end
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
ylim([0.7 1.1])
%if saveFig == 1
%     saveas(fig,'EllipseAreas.png')
% end

%% Find Angles

% this section finds the orientation of the ellipses, as I found that the
% in-built angle from fitellipse.m was inconsistent. the angle found is the
% angle of elevation of the major axis from the positive x axis.

for i = 1:10
    load(['Bivariate_' num2str(i) '_Upper.mat']);
    [a,b] = fitellipse(elpt(:,1),elpt(:,2));
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
    
    load(['Bivariate_' num2str(i) '_Lower.mat']);
    [a,b] = fitellipse(elpt(:,1),elpt(:,2));;
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
ylim([105 140])
if saveFig == 1
    saveas(fig,'EllipseStats.png')
end
