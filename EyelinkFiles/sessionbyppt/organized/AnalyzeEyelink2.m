% Initalize and extract file names
files = dir('*_organized.mat');

pptnames = {};
nmCount = 1;
for i = 1:length(files) 
    if ~any(strcmp(pptnames,extractBefore(files(i).name,'_')))
        pptnames{nmCount} = extractBefore(files(i).name,'_');
        nmCount = nmCount+1;
    end
end

% logical flags for saving figures

saveScatter = true
statsperppt = false
statsforall = true


% Loop through for each participant
for i = 1:length(pptnames)
    clear organized_data
    load([pptnames{i},'_organized.mat']);    
    fig = figure();
    tlo = tiledlayout('flow','TileSpacing','compact');

    % Initialize arrays to store mean distances for upper and lower trials
     upper_mean_dist = zeros(1, length(organized_data.session));
     lower_mean_dist = zeros(1, length(organized_data.session));

    % Loop through for each session 
    for j = 1:length(organized_data.session);

        % Load the upper vs lower hemispheres structure
        load("Hemisphere_Extract.mat");
        Hemispheres = HemispheresByPPT.(pptnames{i}).(append(pptnames{i},'_',num2str(j)));

        % Find average of x and y coordinates from each trial
        x_mean_Lower = [];
        y_mean_Lower = [];
        x_mean_Upper = [];
        y_mean_Upper = [];
        Lower_distance = [];
        Upper_distance = [];
        for k = 1:length(organized_data.session(j).trial);
            if strcmp(Hemispheres{k}, 'Hemisphere: Lower')
                meanx = mean(organized_data.session(j).trial(k).x);
                meany = mean(organized_data.session(j).trial(k).y);
                x_mean_Lower(end+1) = meanx;
                y_mean_Lower(end+1) = meany;
                Lower_distance(end+1) = sqrt((meanx - 720).^2 + (meany - 540).^2);   
            elseif strcmp(Hemispheres{k}, 'Hemisphere: Upper')
                meanx = mean(organized_data.session(j).trial(k).x);
                meany = mean(organized_data.session(j).trial(k).y);
                x_mean_Upper(end+1) = meanx;
                y_mean_Upper(end+1) = meany;
                Upper_distance(end+1) = sqrt((meanx - 720).^2 + (meany - 540).^2); 
            end
        end
        
    
        % scatter plot
        nexttile;
        scatter(x_mean_Lower, y_mean_Lower, 10, [0.9290 0.6940 0.1250]); hold on
        scatter(x_mean_Upper, y_mean_Upper, 10, [0 0.4470 0.7410]); hold on
        scatter(720,540,50,'red', '*', 'LineWidth',2); hold on
        scatter(mean(x_mean_Lower), mean(y_mean_Lower), [0.3 0.3 0.3], '+', 'LineWidth', 1.5);hold on
        scatter(mean(x_mean_Upper), mean(y_mean_Upper), [0.3 0.3 0.3], '+', 'LineWidth', 1.5);hold on
        scatter(480,360,50,[0.3 0.3 0.3]); hold on
        scatter(960,720,50,[0.3 0.3 0.3]); hold on
        scatter(960,360,50,[0.3 0.3 0.3]); hold on
        scatter(480,720,50,[0.3 0.3 0.3]);
        
        %stderror_x_Lower = std(x_mean_Lower) / sqrt(length(x_mean_Lower));
        %stderror_y_Lower = std(y_mean_Lower) / sqrt(length(y_mean_Lower));
        %dx_Lower = stderror_x_Lower * ones(size(x_mean_Lower));
        %dy_Lower = stderror_y_Lower * ones(size(y_mean_Lower));
        %errorbarxy(x_mean_Lower, y_mean_Lower, dx_Lower, dy_Lower);
        
        %stderror_x_Upper = std(x_mean_Upper) / sqrt(length(x_mean_Upper));
        %stderror_y_Upper = std(y_mean_Upper) / sqrt(length(y_mean_Upper));
        %dx_Upper = stderror_x_Upper * ones(size(x_mean_Upper));
        %dy_Upper = stderror_y_Upper * ones(size(y_mean_Upper));
        %errorbarxy(x_mean_Upper, y_mean_Upper, dx_Upper, dy_Upper)
    
      title(sprintf(strcat(['Session',' ',num2str(j)])));
      xlim([0 1440])
      ylim([0 1080])
      xlabel('Screen Width')
      ylabel('Screen Height')

      % Compute mean distance of data points from [720,540] for each session
        lower_mean_dist(j) = nanmean(Lower_distance);
        upper_mean_dist(j) = nanmean(Upper_distance);
        overall_mean_dist(j) = nanmean(cat(2,Lower_distance, Upper_distance));

        lower_mean_y(j) = nanmean(y_mean_Lower);
        upper_mean_y(j) = nanmean(y_mean_Upper);

    end

     sgtitle(strcat(['Average Gaze per Trial:',' ',pptnames{i}]));
     lgd = legend({'Upper','Lower'}, ...
          'FontSize',12);
       legend('boxoff')
       lgd.Layout.Tile = 'east';
       
    if saveScatter
        saveas(gcf,[pwd,'/figures/',strcat(pptnames{i},'_', 'hemisphere','.jpg')]);
    end

    % save stats for participants in an array
    allppt_meand_lower(i) = mean(lower_mean_dist);
    allppt_meand_upper(i) = mean(upper_mean_dist);
    allppt_meand_overall(i) = mean(overall_mean_dist);

    allppt_meany_lower(i) = mean(lower_mean_y);
    allppt_meany_upper(i) = mean(upper_mean_y);

    if statsperppt
        % plot bar chart for each ppt
    
        subplot(1,3,3)
    
        UppervsLower = mwwtest(upper_mean_dist,lower_mean_dist);
    
        b = bar(1:3,[mean(upper_mean_dist), mean(lower_mean_dist), mean(overall_mean_dist)],0.7,'FaceColor',[0.3 0 1]);
        set(gca,'XTickLabel',["Upper","Lower","Overall"]);hold on
        errorbar([mean(upper_mean_dist), mean(lower_mean_dist), mean(overall_mean_dist)],[std(upper_mean_dist)./sqrt(length(upper_mean_dist)), std(lower_mean_dist)./sqrt(length(lower_mean_dist)),std(overall_mean_dist)./sqrt(length(overall_mean_dist))],'LineStyle','none','Color',[0 0 0],'LineWidth',1.5)
        axis square
        ax = gca;
        ax.LineWidth = 2;
        ax.FontName = 'Ariel';
        ax.FontSize = 12;
        title('Mean Distance from Fixation Cross')
        ylabel('Distance (pixels)')
        legend(sprintf('Upper versus Lower, p = %s',string(round(UppervsLower.p(1),4))),'Standard Error')
        ylim([10 70])
    
        saveas(gcf,[pwd,'/figures/',strcat(pptnames{i},'_', 'hemisphereStats','.jpg')]);
    end
end

 if statsforall
        % plot bar chart for all

        subplot(1,3,3)
    
        UppervsLower = mwwtest(allppt_meany_upper,allppt_meany_lower);
    
        b = bar(1:2,[mean(allppt_meany_upper), mean(allppt_meany_lower)],0.7,'FaceColor',[0.3 0 1]);
        set(gca,'XTickLabel',["Upper","Lower","Overall"]);hold on
        errorbar([mean(allppt_meany_upper), mean(allppt_meany_lower)],[std(allppt_meany_upper)./sqrt(length(allppt_meany_upper)), std(allppt_meany_lower)./sqrt(length(allppt_meany_lower))],'LineStyle','none','Color',[0 0 0],'LineWidth',1.5)
        axis square
        ax = gca;
        ax.LineWidth = 2;
        ax.FontName = 'Ariel';
        ax.FontSize = 12;
        title('Mean y from Fixation Cross')
        ylabel('Mean y value')
        ylim([400 600])
        xlim=get(gca,'xlim');
        hold on
        plot(xlim,[540 540],'LineWidth',2,'Color','red')
        legend(sprintf('Upper vs Lower, p = %s',string(round(UppervsLower.p(1),4))),'Standard Error','y value of Fixation Cross')
    
        saveas(gcf,[pwd,'/figures/',strcat('allppt_', 'hemisphereStats','.jpg')]);
    end