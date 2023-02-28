
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
    pptFiles = dir([pptnames{i},'_*.mat']);
    clear session
    if length(pptFiles) > 1
        for y = 1:length(pptFiles)
            load(pptFiles(y).name);
            session(y).gazeloc(:,1) = matFile.Samples.time(:,:);
            session(y).gazeloc(:,2) = matFile.Samples.gx(:,1);
            session(y).gazeloc(:,3) = matFile.Samples.gy(:,1);
            
            for j = 1:length(matFile.Events.Start.time)
                start = matFile.Events.Start.time(:,j);
                stop = matFile.Events.End.time(:,j);
                trial{j} = [];
                trial{j} = [start:stop, trial{j}];
            end
            save([pwd,'/sessionbyppt/' pptnames{i} '_all'],"session")
        end
    end
end
