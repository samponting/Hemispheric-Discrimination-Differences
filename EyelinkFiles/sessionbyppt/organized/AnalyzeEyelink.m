% Initalize and extract file names

files = dir('*.mat');

pptnames = {};
nmCount = 1;
for i = 1:length(files) 
    if ~any(strcmp(pptnames,extractBefore(files(i).name,'_')))
        pptnames{nmCount} = extractBefore(files(i).name,'_');
        nmCount = nmCount+1;
    end
end


for i = 1:length(pptnames)
    clear organized_data
    load([pptnames{i},'_organized.mat']);

    numsessions = length(organized_data.session);
    
    fig = figure()
    tlo = tiledlayout('flow','TileSpacing','compact')
    
    for j = 1:numsessions
    numtrials = length(organized_data.session(j).trial);

    % find average of x and y coordinates from each trial
    x_mean = [];
    y_mean = [];
        for k = 1:numtrials
        meanx = mean(organized_data.session(j).trial(k).x);
        meany = mean(organized_data.session(j).trial(k).y);
        x_mean(k) = meanx;
        y_mean(k) = meany;
        end

    % scatter plot
    nexttile
    scatter(x_mean,y_mean,10,'blue'); hold on
    scatter(720,540,50,'red', '*', 'LineWidth',2); hold on
    scatter(mean(x_mean), mean(y_mean), [0.3 0.3 0.3],'+','LineWidth',1.5);hold on
    scatter(480,360,50,[0.3 0.3 0.3]); hold on
    scatter(960,720,50,[0.3 0.3 0.3]); hold on
    scatter(960,360,50,[0.3 0.3 0.3]); hold on
    scatter(480,720,50,[0.3 0.3 0.3]);
    
    stderror_x = std(x_mean)/sqrt(length(x_mean));
    stderror_y = std(y_mean)/sqrt(length(y_mean));
    dx = stderror_x*ones(size(x_mean));
    dy = stderror_y*ones(size(y_mean));
    errorbarxy(x_mean,y_mean,dx,dy);
    
    title(sprintf(strcat(['Session',' ',num2str(j)])))
    saveas(gcf,[pwd,'/figures/',strcat(pptnames{i},'_','session','_',num2str(j),'.jpg')])

 
    end
end