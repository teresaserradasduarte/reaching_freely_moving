 %% Log files from freely moving
clear; close all; clc

%% % Paths for each mouse and initializaition variables
raw_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\raw_data';
mat_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\analyzed_data\mat_files';
out_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\analyzed_data\output_files';

%data_pathA = '/Users/teresaserradasduarte/Dropbox (Learning Lab)/Learning Lab Team Folder/Patlab protocols/data/TD/raw_data';
project_name = '20211124_AtaxicMice_G3';
%setup = 'freely_mov';
setup = 'freely_moving';
data_path = strcat(raw_folder,filesep,project_name,filesep,setup);
% Find mouse folders
folders_mice = dir(data_path);
folders_mice=folders_mice(4:end,:);
num_animals=size(folders_mice,1);
% Allocate space for mouse cells
mice = cell(num_animals,1);
mice_path = cell(num_animals,1);

% Log file search
searchstr = '*GlobalLogInt';
sessions_TOT = {'S1';'S2'};
nr_sess_all = length(sessions_TOT);

plot_group_total_flag = 1;
if plot_group_total_flag == 1
    ablat = 'dtx';
    ctr_range = [1,2];
    ablat_range = [3,4];
    if strcmp(ablat,"a2a")
        clrs_abl = [216,27,96]./256;
    elseif strcmp(ablat,"dtx")
        clrs_abl = [234,162,33]./256;
    end
    wts_clr = [74,98,116]./256;
end

% Save path
% Save mat file
mat_path=strcat(mat_folder,filesep,project_name,filesep,setup);
mkdir(mat_path);
% Save output files
out_path=strcat(out_folder,filesep,project_name,filesep,setup);
mkdir(out_path);

%% Loop through mice: all sessions
%m=1

% Allocate spare for searching the total nr of trials
nr_trials_all_m = zeros(nr_sess_all,num_animals);


for m = 1: num_animals
    mice(m,1)=cellstr(convertCharsToStrings(folders_mice(m,1).name));
    mice_path(m,1) = strcat(data_path, filesep, mice(m,1));
    
    % Find folders (1 per session)
    folders_session = dir(char(mice_path(m,1)));
    folders_session = folders_session(3:end,:);
    num_sessions = size(folders_session,1);
    % Allocate space for session cells
    sessions = cell(num_sessions,1);
    sessions_path = cell(num_sessions,1);
    
    % Allocate space for log info
    logs = nan(1500,num_sessions);
    timelog =  nan(1500,num_sessions);
    trials_all_ind = nan(400,num_sessions);
    nr_trials_all = zeros(num_sessions,1);
    nr_trials_contigent = zeros(num_sessions,1);
    paw_searching_ind = nan(400,num_sessions);
    trials_contigent_ind = nan(400,num_sessions);
    trials_given_ind = nan(400,num_sessions);
    ITI_duration = nan(400,num_sessions);
    % duration and time
    session_duration = zeros(num_sessions,1);
    duration_reach = nan(400, num_sessions);
    reach_time = nan(400, num_sessions);
    % Reaches per minute
    count_reaches_per_unit_allSess = nan(205, num_sessions);
    units_slide_allSess= nan(205, num_sessions);
    med_units = zeros(num_sessions,1);
    count_reaches_allSess_inRow = [];
    
    %% Load log file for each session
    for s=1:num_sessions
        
        sessions(s,1)=cellstr(convertCharsToStrings(folders_session(s,1).name));
        sessions_path(s,1) = strcat(data_path, filesep, mice(m,1), filesep, sessions(s,1), filesep);
        log_path = wildcardsearch( sessions_path(s,1), searchstr, false,false);
        
        % read log file and write on it
        read_log = [];
        for i=1:size(log_path,1)
            read_log_temp2 = csvread(char(log_path(i,:)));
            read_log_temp = read_log_temp2(2:end,:);
            read_log = cat(1,read_log,read_log_temp);
        end
        logs(1:size(read_log,1),s) = read_log(:,2);
        timelog(1:size(read_log,1),s) = read_log(:,1)-read_log(1,1);
        session_duration(s) = read_log(end,1)-read_log(1,1);
        
        % Interesting parameters
        trials_all_ind_temp = find(logs(:,s)==3);
        trials_all_ind(1:size(trials_all_ind_temp,1),s)=trials_all_ind_temp;
        nr_trials_all(s) = length(trials_all_ind_temp);
        
        paw_searching_ind_temp = find(logs(:,s)==2);
        paw_searching_ind(1:size(paw_searching_ind_temp,1),s)=paw_searching_ind_temp;
        
        % trial type
        if ~isempty(paw_searching_ind_temp)
            is_contigent = 1;
            contigent_or_given_temp = trials_all_ind_temp-paw_searching_ind_temp(1:end-1);
            trials_contigent_ind_temp = trials_all_ind_temp(contigent_or_given_temp==1);
            trials_contigent_ind(1:size(trials_contigent_ind_temp,1),s)=trials_contigent_ind_temp;
            trials_given_ind_temp = trials_all_ind_temp(contigent_or_given_temp>1);
            trials_given_ind(1:size(trials_given_ind_temp,1),s)=trials_given_ind_temp;
            nr_trials_contigent(s) = length(trials_contigent_ind_temp);
        else
            is_contigent = 0;
            water_reached =find(logs(:,s)==31); 
        end
        
        % ITIs duration
        ITIs_stop_ind=find(logs(:,s)==40);
        ITI_duration_temp = timelog(ITIs_stop_ind,s)-timelog(ITIs_stop_ind-1,s);
        ITI_duration(1:size(ITI_duration_temp,1),s)=ITI_duration_temp;
        
        
        % Reach duration ( = reaction time + movement time + contigency time)
        if is_contigent==1
        duration_reach(1:size(trials_contigent_ind_temp,1),s) = timelog(trials_contigent_ind_temp,s)-timelog(trials_contigent_ind_temp-1,s);
        % Time of reach
        reach_time(1:size(trials_contigent_ind_temp,1),s) = timelog(trials_contigent_ind_temp,s);
        else
            duration_reach(1:size(water_reached,1),s) = timelog(water_reached,s)-timelog(water_reached-1,s);
            reach_time(1:size(water_reached,1),s) = timelog(water_reached,s);
        end

        
        % Counts per unit time (sliding window through session)
        sliding_win = 200;
        sliding_step = 10;
        max_ind = floor((session_duration(s)-sliding_win)/sliding_step);
        count_reaches_per_unit = [];
        units_slide = [];
        for i = 1:max_ind
            sliding_vec = i*sliding_step+1:i*sliding_step+sliding_win;
            reaches_current_unit = numel(find(reach_time(:,s)>sliding_vec(1) & reach_time(:,s)<sliding_vec(end)));
            units_slide = [units_slide;sliding_vec(end)];
            count_reaches_per_unit=[count_reaches_per_unit;reaches_current_unit];
        end
        count_reaches_per_unit_allSess(1:size(count_reaches_per_unit,1),s)=count_reaches_per_unit;
        units_slide_allSess(1:size(units_slide,1),s) = units_slide;
        med_units(s)=round(numel(units_slide)/2);
        % Count in row
        count_reaches_allSess_inRow = [count_reaches_allSess_inRow;nan(10,1);count_reaches_per_unit];
        time_inRow = linspace(0,numel(count_reaches_allSess_inRow)*(sliding_step/60),numel(count_reaches_allSess_inRow));


        
        % PLOTS --------------------------------------
        save_sessDir=char(strcat(out_path,filesep,mice(m,1),filesep,sessions(s,1)));
        mkdir(save_sessDir)
        %
        % Plot durations: histogram and durations over time
        figure(1)
        subplot(121)
        nbins=200;
        %     maxval = max(duration_reach(:,si)); minval = min(duration_reach(:,si));
        %     binstep = (maxval-minval)/(nbins-1);
        %     bins = minval:binstep:maxval; % from minval to maxval with binstep increments (xx scale)
        %     histogram(duration_reach(:,si),bins,'normalization','probability');
        histogram(duration_reach(:,s),'normalization','probability','binwidth',0.5);
        xlim([-1 20])
        ylabel('counts (/total)'); xlabel('reach duration (sec)'); axis square; shg
        % Duration over time
        subplot(122)
        plot(squeeze(reach_time(:,s)),squeeze(duration_reach(:,s)),'.');
        xlabel('time (sec)'), ylabel('duration of reach (sec)'); axis square
        saveas(gcf,strcat(save_sessDir,filesep,'reach_duration.png'),'png')
        pause(1)

        figure(2)
        set(gcf, 'Position', [2069 421 1605 399])
        subplot(121)
        reaches = nan(size(reach_time(:,s))); reaches(~isnan(reach_time(:,s)))=1;
        plot(reach_time(:,s),reaches,'.'); 
        xlim([0 session_duration(s)])
        xlabel('time (sec)'), ylabel('reaches')
         
        subplot(122)
        plot(units_slide,count_reaches_per_unit,'o')
        title(sprintf('%s%d%s','Number of reaches per ',sliding_win,' seconds'))
        ylabel('number of reaches'); xlabel('time (sec)');
        saveas(gcf,strcat(save_sessDir,filesep,'trialsOverTime_slidingWin_reachesCount.png'),'png')    
        pause(1)
    end
    close all
    
    
    %% Save mat file
    mouseX = mice(m,:);
    saveDir_mouseX = char(strcat(mat_path,filesep,mice(m,1)));
    mkdir(saveDir_mouseX)
    save(strcat(saveDir_mouseX,filesep,'logdata.mat'),...
        'mouseX','sessions','sessions_path',...
        'logs','timelog','trials_all_ind','nr_trials_all',...
        'trials_contigent_ind','nr_trials_contigent','trials_given_ind',...
        'ITI_duration','duration_reach','reach_time','session_duration',...
        'sliding_win','sliding_step','count_reaches_per_unit_allSess','units_slide_allSess');
    
    
    %% Plot data
    % distribution of durations
    saveDirOut_mouseX = char(strcat(out_path,filesep,mice(m,1)));
    figure(1)
    subplot(221)
    %nbins=200;
    for ss=1:num_sessions
        %maxval = max(duration_reach(:,ss)); minval = min(duration_reach(:,ss));
        %binstep = (maxval-minval)/(nbins-1);
        %bins = minval:binstep:maxval; % from minval to maxval with binstep increments (xx scale)
        h = histogram(duration_reach(:,ss),'normalization','probability','binwidth',0.2);
        hold on
    end
    hold off
    xlim([-0.2 6])
    ylabel('counts (/total)'); xlabel('reach duration (sec)'); axis square; shg
    legend(sessions)
    colorOrder = get(gca, 'ColorOrder');
    
    subplot(222)
    si=1;
    for ss=1:num_sessions
        h = histogram(duration_reach(:,ss),'normalization','probability','binwidth',0.2,'Visible',0);
        hc = contourhist(h, 'Color', colorOrder(si,:), 'LineWidth', 1.5);
        hold on
        si=si+1;
    end
    hold off
    xlim([-0.2 6])
    ylabel('counts (/total)'); xlabel('reach duration (sec)'); axis square; shg
    
    % Duration over time
    subplot(223)
    plot(squeeze(reach_time),squeeze(duration_reach),'.','markersize',15);
    xlabel('time (sec)'), ylabel('duration of reach (sec)'); axis square;
    legend(sessions)
    ylim([0 50])
    
    % Mean duration over session
    subplot(224)
    plot(nanmean(duration_reach),'-ok')
    xlabel('session'), ylabel('mean duration'); axis square;
    axis([0.9 num_sessions+.1 0 10])
    saveas(gcf,strcat(saveDirOut_mouseX,filesep,'duration_reach.png'),'png')
    
    
    % ------------------------------
    % Number of trials across sessions
    figure(2)
    set(gcf, 'Position', [2069 421 1605 399])
    plot(time_inRow,count_reaches_allSess_inRow,'-ko');
    set(gca,'linewidth',1.5,'box','off','layer','top','GridAlpha',0.05,'TickLength',[0.002, 0.001]);
    xline(time_inRow(diff(isnan(count_reaches_allSess_inRow))==1),'--','LineWidth',2,'Color',[.8 .8 .8])

    hold on
    %plot(tickss(2:4:end),nr_trials_all/10,'r*')
    hold off
    ylim([0 50])
    title(sprintf('%s%d%s','Number of trials per ',sliding_win,' seconds, across sessions'))
    
    %     xticks(tickss);
    %     xticklabels(tick_labels)
    ylabel('Number of water droplets delivered')
    xlabel('time across sessions (min)')
    saveas(gcf,strcat(saveDirOut_mouseX,filesep,'number_trials_overSessions_slidingWin.png'),'png')
    saveas(gcf,strcat(saveDirOut_mouseX,filesep,'number_trials_overSessions_slidingWin.epsc'),'epsc')


    %%
    figure
    plot(nr_trials_all,'-ko','LineWidth',1.5);
    set(gca,'linewidth',1.5,'box','off','layer','top','GridAlpha',0.05,'TickLength',[0.002, 0.001]);
    title('Nr of trials per session')
    xticklabels(sessions)
    xlabel('sessions'); ylabel('number of trials')
    saveas(gcf,strcat(saveDirOut_mouseX,filesep,'number_trials_session.png'),'png')


    nr_trials_all_m(:,m) = nr_trials_all;

end


%%
 nxh1 = 1;
 max_n_group = max(length(ctr_range),length(ablat_range));
 vec_nr_trials =nan(max_n_group,(nxh1*3)-1);
 vec_nr_trials(1:length(ctr_range),1:3:nxh1*3)=nr_trials_all_m(2,ctr_range)';
 vec_nr_trials(1:length(ablat_range),2:3:nxh1*3+1)=nr_trials_all_m(2,ablat_range)';
 length_n_Trials=size(vec_nr_trials,1);
 [h_nTrials,sig_nTrials]=ttest2(vec_nr_trials(:,1),vec_nr_trials(:,2));


 %%

 clrs = cat(1,wts_clr,clrs_abl);
 figure()

boxplot(vec_nr_trials,'Symbol', 'k.','Color','k','Widths',0.8);
hold on
for pos=1:size(vec_nr_trials,2)
    if ismember(pos,1:3:nxh1*3), c=1; else, c=2; end
    f=scatter(ones(1,length_n_Trials)+(pos-1).*(1+(rand(01,length_n_Trials)-0.5)/50),vec_nr_trials(:,pos),[],clrs(c,:),'filled','LineWidth',1.5);
end
hold off
axis square;% ylim([-20 250])
XTickLabel={'CTR';'DTX'};
set(gca,'XTickLabel',XTickLabel); ylabel('Nr of trials');
set(gca,'linewidth',1.2,'box','off','layer','top','GridAlpha',0.05);
title({'Nr trials',sprintf('%s%1.3f','sig = ',sig_nTrials)})

saveas(gcf,strcat(out_path,filesep,'boxplot_nr_trials.png'),'png')
saveas(gcf,strcat(out_path,filesep,'boxplot_nr_trials.epsc'),'epsc')
    
    


 
%%

axeOpt = {'linewidth',1.5,'box','off','GridAlpha',0.05,'ticklength',[1,1]*.01};
lw = 1.5;

figure()
plot(nr_trials_all_m(:,ctr_range),'LineWidth',lw,'Color',wts_clr); hold on
plot(nr_trials_all_m(:,ablat_range),'LineWidth',lw,'Color',clrs_abl);hold off
xlim([.8 nr_sess_all+.2])
set(gca,axeOpt{:})
title('Nr of trials per session')
xticks(1:nr_sess_all)
xticklabels(sessions_TOT)
xlabel('sessions'); ylabel('number of trials')
legend(cat(1,mice(ctr_range),mice(ablat_range)),'box','off','location','southeastoutside','Interpreter', 'none')
saveas(gcf,strcat(out_path,filesep,'number_trials_session_learning.png'),'png')
saveas(gcf,strcat(out_path,filesep,'number_trials_session_learning.epsc'),'epsc')
    
    

%% Save
save(strcat(mat_path,filesep,'logdata_load.mat'),...
    'mice','folders_mice','sessions','folders_session',...
    'data_path','mat_path','out_path','sessions','nr_trials_all_m',...
    'vec_nr_trials','sig_nTrials');
%
%     'duration_reach','reach_time','tickss','tick_labels',...
%     'count_reaches_per_unit_allSess','sliding_win','nr_trials_all');

shg







