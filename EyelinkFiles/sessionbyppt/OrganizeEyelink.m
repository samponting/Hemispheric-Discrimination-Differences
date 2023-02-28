% names = {'SJP' 'MLV' 'LW' 'LR' 'JLA' 'JG' 'DJT'};

files = dir('*.mat');

pptnames = {};
nmCount = 1;
for i = 1:length(files) 
    if ~any(strcmp(pptnames,extractBefore(files(i).name,'_')))
        pptnames{nmCount} = extractBefore(files(i).name,'_');
        nmCount = nmCount+1;
    end
end

for lol = 1:length(pptnames)
    
    load([pptnames{lol},'_all.mat'])

    organized_data = struct('session', struct('trial', struct('time', {}, 'x', {}, 'y', {})));
    
    for s = 1:length(session)
        
        % Helper variables to separate out each trial
        organized_data.session(s) = struct('trial', struct('time', {}, 'x', {}, 'y', {}));
        session_data = session(s).gazeloc;
        prev_time = session_data(1, 1);
        times = [];
        xs = [];
        ys = [];
        
        % For each trial
        for i = 1:length(session_data)
            curr_time = session_data(i, 1);
            if curr_time - prev_time > 1
                % A jump in time occurred, indicating the start of a new recording
                organized_data.session(s).trial(end + 1) = struct('time', times, 'x', xs, 'y', ys);
                times = [];
                xs = [];
                ys = [];
            else
                % Record trial data row temporarily to be later added at once
                times(end + 1) = session_data(i, 1);
                xs(end + 1) = session_data(i, 2);
                ys(end + 1) = session_data(i, 3);
            end
            prev_time = curr_time;
        end
    
    end

    save([pwd,'/organized/' pptnames{lol} '_organized'],"organized_data")

end