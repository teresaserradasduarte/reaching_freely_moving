%% Log files from freely moving
clear; close all; clc

%% % Paths for each mouse and initializaition variables
raw_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\raw_data';
mat_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\analyzed_data\mat_files';
out_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\analyzed_data\output_files';

%data_pathA = '/Users/teresaserradasduarte/Dropbox (Learning Lab)/Learning Lab Team Folder/Patlab protocols/data/TD/raw_data';
project_name = '20250106_A2aCasp_G2';
setup = 'freely_mov';
data_path = strcat(mat_folder,filesep,project_name,filesep,setup);
save_path = strcat(out_folder,filesep,project_name,filesep,setup,filesep,'group_results',filesep,'general');
if ~exist('save_path','dir'), mkdir(save_path); end

% Find mouse folders
folders_mice = dir(data_path);
folders_mice=folders_mice(3:end,:);
num_animals=size(folders_mice,1);
% Allocate space for mouse cells
mice = cell(num_animals,1);
mice_path = cell(num_animals,1);
num_sessions_tmp = 2;

% Allocate space for variables
fps = 120;
sess_dur_long = fps*60*45;
nr_datapoints_trials = fps*60*6;
nr_trials_max = 400;
nr_datapoints = 101;
nr_reaches_max = 1600;

paw_ss_all = nan(sess_dur_long,2,num_sessions_tmp,num_animals);
paw_tt_all =  nan(nr_datapoints_trials,2,nr_trials_max,num_sessions_tmp,num_animals);
snout_tt_all = nan(nr_datapoints_trials,2,nr_trials_max,num_sessions_tmp,num_animals);
tongue_tt_all = nan(nr_datapoints_trials,2,nr_trials_max,num_sessions_tmp,num_animals);
reaches_rr_all = nan(nr_datapoints,2,nr_reaches_max,num_sessions_tmp,num_animals);


for m = 1: num_animals
    mice(m,1)=cellstr(convertCharsToStrings(folders_mice(m,1).name));
    mice_path(m,1) = strcat(data_path, filesep, mice(m,1));
    % Find folders (1 per session)
    folders_session = dir(char(mice_path(m,1)));
    folders_session = folders_session(5:end-1,:); % RANGE ADAPTD TO SESISONS OF INTEREST
    num_sessions = size(folders_session,1);
    % Allocate space for session cells
    sessions_name = cell(num_sessions,1);
    sessions_path = cell(num_sessions,1);

    for s=1:num_sessions
        sessions_name(s,1)=cellstr(convertCharsToStrings(folders_session(s,1).name));
        sessions_path(s,1) = strcat(data_path, filesep, mice(m,1), filesep, sessions_name(s,1), filesep);

        mat_file1 = strcat(char(sessions_path(s,1)),'session_reaching_data.mat');
        if exist(mat_file1,"file")
            fprintf('%s%s%s%s%s','loading mouse ',char(mice(m,1)),', session ',char(sessions_name(s,1)))
            fprintf('\n')

            load(mat_file1,'trials','sessions','reaches')

            paw_ss_all(1:size(sessions.paw_ss,1),:,s,m) = sessions.paw_ss;
            paw_tt_all(1:size(trials.paw,1),:,1:size(trials.paw,3),s,m) = trials.paw;
            snout_tt_all(1:size(trials.paw,1),:,1:size(trials.paw,3),s,m) = trials.snout;
            tongue_tt_all(1:size(trials.paw,1),:,1:size(trials.paw,3),s,m) = trials.tongue;
            reaches_rr_all(:,:,1:size(reaches.reach_mat,3),s,m) = reaches.reach_mat;
        else
            fprintf('%s%s%s%s','no data for mouse ',char(mice(m,1)),', session ',char(sessions_name(s,1)))
            fprintf('\n')

        end
    end
end


%% Plot stuff

a2a_range = [];
wts_range = [];
for m = 1:num_animals
    phenotype = a2a_phenotype(mice{m});
    if strcmp(phenotype,'ctr')
        wts_range = [wts_range,m];
    elseif strcmp(phenotype,'a2a')
        a2a_range = [a2a_range,m];
    end
end

a2a_clr = [216,27,96]./256;
wts_clr = [74,98,116]./256;
clrs_m = repmat(wts_clr,[num_animals,1]);
clrs_m(a2a_range,:) = repmat(a2a_clr, [length(a2a_range),1]);



cat1_clr=[36 62 54]./256;
cat2_clr=[124 169 130]./256;
cat3_clr=[224 238 198]./256;
%indiv_a_clr = [.8 .8 .8];


axeOpt = {'linewidth',1.5,'box','off','GridAlpha',0.05,'ticklength',[1,1]*.01};


%% Check all sessiom

time_sess = (0:size(paw_ss_all,1)-1)./fps;
time_selec =  3.5*60*fps:4*60*fps;
x_lim = [80 350];
z_lim = [100 350];
pos_rec_x = [time_sess(time_selec(1))./60 x_lim(1) ...
    (time_sess(time_selec(end))./60)-(time_sess(time_selec(1))./60)...
    x_lim(end)-x_lim(1)];
pos_rec_z = [time_sess(time_selec(1))./60 z_lim(1) ...
    (time_sess(time_selec(end))./60)-(time_sess(time_selec(1))./60)...
    z_lim(end)-z_lim(1)];

s=3;
for m=1:num_animals
    figure();
    tt = tiledlayout(4,1);

    title(tt,sprintf('%s%s%s',char(mice(m)), ' , ',char(sessions_name(s))),'Interpreter','none')

    nexttile
    plot(time_sess./60,paw_ss_all(:,1,s,m),'color',clrs_m(m,:));
    rectangle('Position',pos_rec_x, 'EdgeColor',[.8 .8 .8],'LineWidth',2)
    set(gca,axeOpt{:})
    ylim([x_lim])
    xlim([0 30])
    xlabel('time in session (min)');
    ylabel('x (px)');

    nexttile
    plot(time_sess./60,paw_ss_all(:,2,s,m),'color',clrs_m(m,:));
    rectangle('Position',pos_rec_z, 'EdgeColor',[.8 .8 .8],'LineWidth',2)
    set(gca,axeOpt{:})
    ylim([z_lim])
    xlim([0 30])
    xlabel('time in session (min)');
    ylabel('z (px)');


    nexttile
    plot(time_sess(time_selec)./60,paw_ss_all(time_selec,1,s,m),'color',clrs_m(m,:));
    set(gca,axeOpt{:})
    ylim([x_lim])
    xlim([time_sess(time_selec(1))./60 time_sess(time_selec(end))./60])
    xlabel('time (min)');
    ylabel('x (px)');

    nexttile
    plot(time_sess(time_selec)./60,paw_ss_all(time_selec,2,s,m),'color',clrs_m(m,:));
    set(gca,axeOpt{:})
    ylim([x_lim])
    xlim([time_sess(time_selec(1))./60 time_sess(time_selec(end))./60])
    xlabel('time  (min)');
    ylabel('z (px)');

    set(gcf,'position',[1925         123        1558         863],'Color','w')

    saveas(gcf,strcat(save_path,filesep,char(mice(m)),'_',char(sessions_name(s)),'_pawSessT.png'),'png')

end


%%
%reaches_rr_all(reaches_rr_all==0)=nan;
transpa = 0.03;
x_lim = [0 350];
z_lim = [0 450];
%s=4;
for s=1:num_sessions
    figure();
    tt = tiledlayout(2,num_animals/2);
    title(tt,sprintf('%s%s','reach trajectories ',char(sessions_name(s))),'Interpreter','none')
    for m = 1:num_animals
        nexttile
        plot(squeeze(reaches_rr_all(:,1,:,s,m)),...
            squeeze(reaches_rr_all(:,2,:,s,m)),"Color",cat(2,clrs_m(m,:),transpa));
        hold on
        %          plot(median(reaches_rr_all(:,1,:,s,m),3,'omitnan'),...
        %              median(reaches_rr_all(:,2,:,s,m),3,'omitnan'),...
        %              "Color",clrs_m(m,:),'LineWidth',2);
        set(gca,axeOpt{:})
        axis square;
        xlabel('x'); ylabel('y')
        axis([x_lim z_lim])
        title(char(mice(m)),'Interpreter','none')
    end
    set(gcf,'Color','w', 'position',[1975         161        1812         821])
    saveas(gcf,strcat(save_path,filesep,'reach_trajectories_',char(sessions_name(s))),'png');
end


%% Map of reaches

s=4;
nbins = 50;
for s=1:num_sessions
    figure();
    tt = tiledlayout(2,num_animals/2);
    title(tt,sprintf('%s%s','reach trajectories ',char(sessions_name(s))))
    for m = 1:num_animals
        nexttile

    histogram2(squeeze(reaches_rr_all(:,1,:,s,m)),...
        squeeze(reaches_rr_all(:,2,:,s,m)),...
        nbins,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none',...
        'normalization','pdf');
    set(gca,'CLim',[0 7E-5])
    colorbar
    %set(gca, 'XTick', [], 'YTick', [],'box','off')
    grid off
        axis square;
        xlabel('x'); ylabel('y')
        axis([0 300 0 450])
        title(char(mice(m)),'color',clrs_m(m,:),'Interpreter','none')
    end
    grid off
    set(gcf,'Color','w', 'position',[1975         161        1812         821])
    saveas(gcf,strcat(save_path,filesep,'reach_occupancy_',char(sessions_name(s))),'png');
end


%% Paw all occupancy
nbins = 50;
for s=1:num_sessions
    figure();
    tt = tiledlayout(2,num_animals/2);
    title(tt,sprintf('%s%s','reach trajectories ',char(sessions_name(s))))
    for m = 1:num_animals
        nexttile

    histogram2(squeeze(paw_ss_all(:,1,s,m)),...
        squeeze(paw_ss_all(:,2,s,m)),...
        nbins,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none',...
        'normalization','pdf');
    set(gca,'CLim',[0 2E-5])
    colorbar
    %set(gca, 'XTick', [], 'YTick', [],'box','off')
        axis square;
        xlabel('x'); ylabel('y')
        axis([0 300 0 450])
        title(char(mice(m)),'color',clrs_m(m,:),'Interpreter','none')
            grid off
    end
    set(gcf,'Color','w', 'position',[1975         161        1812         821])
  
    saveas(gcf,strcat(save_path,filesep,'paw_all_session_occupancy',char(sessions_name(s))),'png');
end

%% DIFFERENCE PAW_ALL_OCCUPACY AND REACHES OCCUPANCY 

%% Snout all occupancy
%s=4;
nbins = 50;
for s=1:num_sessions
    figure();
    tt = tiledlayout(2,num_animals/2);
    title(tt,sprintf('%s%s','reach trajectories ',char(sessions_name(s))))
    for m = 1:num_animals
        nexttile

    histogram2(squeeze(snout_tt_all(:,1,:,s,m)),...
        squeeze(snout_tt_all(:,2,:,s,m)),...
        nbins,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none',...
        'normalization','pdf');
    set(gca,'CLim',[0 9E-7])
    colorbar
    %set(gca, 'XTick', [], 'YTick', [],'box','off')
        axis square;
        xlabel('x'); ylabel('y')
        %axis([30 300 50 350])
        title(char(mice(m)),'color',clrs_m(m,:),'Interpreter','none')
        grid off
    end
    set(gcf,'Color','w', 'position',[1975         161        1812         821])
    saveas(gcf,strcat(save_path,filesep,'snout_occupancy',char(sessions_name(s))),'png');
end


%% Tongue all occupancy
%s=4;
nbins = 50;
for s=1:num_sessions
    figure();
    tt = tiledlayout(2,num_animals/2);
    title(tt,sprintf('%s%s','tongue occupancy ',char(sessions_name(s))))
    for m = a2a_range

        nexttile
        histogram2(squeeze(tongue_tt_all(:,1,:,s,m)),...
            squeeze(tongue_tt_all(:,2,:,s,m)),...
            nbins,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none',...
            'normalization','pdf');
        set(gca,'CLim',[0 9E-8])
        colorbar
        title(char(mice(m)),'color',clrs_m(m,:),'Interpreter','none')
        %set(gca, 'XTick', [], 'YTick', [],'box','off')
        axis square;
        xlabel('x'); ylabel('y')
        title(char(mice(m)),'color',clrs_m(m,:),'Interpreter','none')
    end
    for m = wts_range
        nexttile
        histogram2(squeeze(tongue_tt_all(:,1,:,s,m)),...
            squeeze(tongue_tt_all(:,2,:,s,m)),...
            nbins,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none',...
            'normalization','pdf');
        set(gca,'CLim',[0 9E-8])
        colorbar
        title(char(mice(m)),'color',clrs_m(m,:),'Interpreter','none')
        %set(gca, 'XTick', [], 'YTick', [],'box','off')
        axis square;
        xlabel('x'); ylabel('y')
        %axis([30 300 50 350])
    end
    set(gcf,'Color','w', 'position',[1975         161        1812         821])
    saveas(gcf,strcat(save_path,filesep,'tongue_occupancy',char(sessions_name(s))),'png');
end
%%
close all



