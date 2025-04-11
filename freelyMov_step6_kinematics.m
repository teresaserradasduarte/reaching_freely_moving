%% Step 6 - kinematics
clear; close all; clc

%% % Paths for each mouse and initializaition variables
raw_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\raw_data';
mat_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\analyzed_data\mat_files';
out_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\analyzed_data\output_files';

%data_pathA = '/Users/teresaserradasduarte/Dropbox (Learning Lab)/Learning Lab Team Folder/Patlab protocols/data/TD/raw_data';
project_name = '20250106_A2aCasp_G2';
setup = 'freely_mov';
data_path = strcat(mat_folder,filesep,project_name,filesep,setup);
save_path = strcat(out_folder,filesep,project_name,filesep,setup,filesep,'group_results');
if ~exist('save_path','dir'), mkdir(save_path); end

% Find mouse folders
folders_mice = dir(data_path);
folders_mice=folders_mice(3:end,:);
num_animals=size(folders_mice,1);
% Allocate space for mouse cells
mice = cell(num_animals,1);
mice_path = cell(num_animals,1);


% Allocate pace for variables
n_sess_TOT = 5;
n_sess = 5;
nr_reaches_max = 1600;
nr_reaches_max_TOT = 1600;
%nr_trials_max = 400;
nr_datapoints = 101;

nr_trials_TOT = nan(n_sess_TOT,num_animals);
nr_trials_t = nan(n_sess,num_animals);
nr_reaches_t = nan(n_sess,num_animals);
length_sess = nan(n_sess,num_animals);
reach_times_sess = nan(nr_reaches_max_TOT,n_sess,num_animals);
dur_reaches_t = nan(nr_reaches_max,n_sess,num_animals);
r_paw_used_t = nan(nr_reaches_max,n_sess,num_animals);
start_fwd_t = nan(nr_reaches_max,n_sess,num_animals);
stop_reach_t = nan(nr_reaches_max,n_sess,num_animals);

cat_reach_id_t = nan(nr_reaches_max,n_sess,num_animals);
class_reaches_t = nan(3,nr_reaches_max,n_sess,num_animals);

iri_sec_t = nan(nr_reaches_max_TOT,n_sess,num_animals);
time_to_1st_reach_peak_t = nan(nr_reaches_max,n_sess,num_animals);
time_to_1st_reach_start_t = nan(nr_reaches_max,n_sess,num_animals);
reaches_all_t = nan(nr_datapoints,2,nr_reaches_max,n_sess,num_animals);
reaches_start_stop_t = nan(nr_datapoints,2,nr_reaches_max,n_sess,num_animals);
speed_reaches_all_t = nan(nr_datapoints,2,nr_reaches_max,n_sess,num_animals);
distance_water_all_t = nan(nr_datapoints,nr_reaches_max,n_sess,num_animals);
speed_water_all_t = nan(nr_datapoints,nr_reaches_max,n_sess,num_animals);

cat_reach_nr_t = nan(3,n_sess,num_animals);
puorpose_yesNo_nr = nan(2,n_sess,num_animals);
hit_yesNo_nr = nan(2,n_sess,num_animals);
suc_yesNo_nr = nan(2,n_sess,num_animals);

for m = 1: num_animals
    mice(m,1)=cellstr(convertCharsToStrings(folders_mice(m,1).name));
    mice_path(m,1) = strcat(data_path, filesep, mice(m,1));
    % Find folders (1 per session)
    folders_session = dir(char(mice_path(m,1)));
    folders_session = folders_session(3:end-1,:); % RANGE ADAPTD TO SESISONS OF INTEREST
    num_sessions = size(folders_session,1);
    % Allocate space for session cells
    sessions_name = cell(num_sessions,1);
    sessions_path = cell(num_sessions,1);

    for s=1:num_sessions
        sessions_name(s,1)=cellstr(convertCharsToStrings(folders_session(s,1).name));
        sessions_path(s,1) = strcat(data_path, filesep, mice(m,1), filesep, sessions_name(s,1), filesep);

        mat_file1 = strcat(char(sessions_path(s,1)),'session_reaching_data.mat');
        mat_file2 = strcat(char(sessions_path(s,1)),'reaches_analysis.mat');

        if exist(mat_file1,"file")
            fprintf('%s%s%s%s%s','loading mouse ',char(mice(m,1)),', session ',char(sessions_name(s,1)))
            fprintf('\n')

            load(mat_file1,'trials','reaches','sessions')
            length_sess(s,m) = size(sessions.paw_ss,1);
            reach_times_sess(1:length(sessions.pks_frames),s,m) = sessions.pks_frames;
            nr_trials_t(s,m) = trials.n_trials;
            nr_reaches_t(s,m) = reaches.r_num;
            dur_reaches_t(1:nr_reaches_t(s,m),s,m) = reaches.start_forw.dur_forw_mat;
            r_paw_used_t(1:nr_reaches_t(s,m),s,m) = reaches.Rpaw_used;
            start_fwd_t(1:nr_reaches_t(s,m),s,m) = reaches.start_forw.start_forw_mat;
            stop_reach_t(1:nr_reaches_t(s,m),s,m) = reaches.start_stop.stop_reach_mat;

            cat_reach_id_t(1:nr_reaches_t(s,m),s,m) = reaches.cat_reach;
            class_reaches_t(1,1:nr_reaches_t(s,m),s,m) = reaches.purpose_reach;
            class_reaches_t(2,1:nr_reaches_t(s,m),s,m) = reaches.hit_reach;
            class_reaches_t(3,1:nr_reaches_t(s,m),s,m) = reaches.success_reach;
            tm = reaches.time_range;
            fps = sessions.video.frame_rate;
            clear trials reaches sessions;

            load(mat_file2,...
                'iri_sec','time_to_1st_reach_peak','time_to_1st_reach_start',...
                'reaches_all','reaches_start_stop','speed_reaches_all',...
                'reach_to_water','speed_reach_to_water',...
                'perc1','perc2','perc3','porpuse_nr','noporpuse_nr',...
                'hit_nr','miss_nr','succ_nr','fail_nr');

            iri_sec_t(1:length(iri_sec),s,m) = iri_sec;

            time_to_1st_reach_peak_t(1:length(time_to_1st_reach_peak),s,m) = time_to_1st_reach_peak;
            time_to_1st_reach_start_t(1:length(time_to_1st_reach_start),s,m) = time_to_1st_reach_start;
            reaches_all_t(:,:,1:nr_reaches_t(s,m),s,m) = reaches_all;
            reaches_start_stop_t(:,:,1:nr_reaches_t(s,m),s,m) = reaches_start_stop;
            speed_reaches_all_t(:,:,1:nr_reaches_t(s,m),s,m) = speed_reaches_all;
            distance_water_all_t(:,1:nr_reaches_t(s,m),s,m) = reach_to_water;
            speed_water_all_t(:,1:nr_reaches_t(s,m),s,m) = speed_reach_to_water;

            cat_reach_nr_t(1,s,m) = perc1;
            cat_reach_nr_t(2,s,m) = perc2;
            cat_reach_nr_t(3,s,m) = perc3;

            puorpose_yesNo_nr(1,s,m) = porpuse_nr;
            puorpose_yesNo_nr(2,s,m) = noporpuse_nr;
            hit_yesNo_nr(1,s,m) = hit_nr;
            hit_yesNo_nr(2,s,m) = miss_nr;
            suc_yesNo_nr(1,s,m) = succ_nr;
            suc_yesNo_nr(2,s,m) = fail_nr;

            clear iri_sec reaches_all speed_reaches_all ...
                reach_to_water speed_reach_to_water ...
                perc1 perc2 perc3 porpuse_nr noporpuse_nr ...
                hit_nr miss_nr succ_nr fail_nr

        else
            fprintf('%s%s%s%s','no data for mouse ',char(mice(m,1)),', session ',char(sessions_name(s,1)))
            fprintf('\n')

        end
    end
    % Load nr_trials of all sessions
    load(char(strcat(data_path, filesep, mice(m,1), filesep,'logdata.mat')),'nr_trials_all');
    nr_trials_TOT(:,m) = nr_trials_all;
    clear nr_trials_all
end
fprintf(['%s' ...
    '\n'],'done!')

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


axeOpt = {'linewidth',1.5,'box','off','GridAlpha',0.05,'ticklength',[1,1]*.01,'TickDir','out'};


%% Speed forward only
reaches_forward_part = nan(41,2,nr_reaches_max,n_sess,num_animals);

for m=1:num_animals
    for s=1:n_sess
        for rr = 1:nr_reaches_max
            if ~isnan(start_fwd_t(rr,s,m))
                reaches_forward_part(start_fwd_t(rr,s,m):41,:,rr,s,m) = reaches_all_t(start_fwd_t(rr,s,m):41,:,rr,s,m);
            else
                reaches_forward_part(:,:,rr,s,m) = nan(41,2);
            end

        end
    end
end
reaches_forward_part(reaches_forward_part==0)=nan;
reaches_fwd_last = squeeze(reaches_forward_part(:,:,:,5,:));
dt = 1/fps;
size_npt = size(reaches_forward_part,1);

%% Instantaneous speed
vel_fwd_partial = diff(reaches_fwd_last,1);
speed_fwd = squeeze(sqrt(vel_fwd_partial(:,1,:,:).^2+vel_fwd_partial(:,2,:,:).^2)./dt);

%% Instantaneous acelleration
acc_fwd_partial = nan(size_npt-2,2,nr_reaches_max,num_animals);
for i=1:size_npt-2
    acc_fwd_partial(i,:,:,:) = reaches_fwd_last(i+2,:,:,:) - 2.*reaches_fwd_last(i+1,:,:,:) + reaches_fwd_last(i,:,:,:);
end
acc_fwd = squeeze(sqrt(acc_fwd_partial(:,1,:,:).^2+acc_fwd_partial(:,2,:,:).^2)./dt.^2);

%% Instantaneous jerk
jerk_fwd_partial = nan(size_npt-3,2,nr_reaches_max,num_animals);
for i=1:size_npt-3
    jerk_fwd_partial(i,:,:,:) = reaches_fwd_last(i+3,:,:,:) - 3*reaches_fwd_last(i+2,:,:,:) ...
        + 3.*reaches_fwd_last(i+1,:,:,:) - reaches_fwd_last(i,:,:,:);
end
jerk_fwd = squeeze(sqrt(jerk_fwd_partial(:,1,:,:).^2+jerk_fwd_partial(:,2,:,:).^2)./dt.^3);

%% Mean vector 

mean_vel = squeeze(mean(mean(speed_fwd,1,'omitnan'),2,'omitnan'));
mean_acc = squeeze(mean(mean(acc_fwd,1,'omitnan'),2,'omitnan'));
mean_jrk = squeeze(mean(mean(jerk_fwd,1,'omitnan'),2,'omitnan'));
nxh3 = 3;
vec_mean_vel = cat(2,mean_vel(wts_range),mean_vel(a2a_range));
vec_mean_acc = cat(2,mean_acc(wts_range),mean_acc(a2a_range));
vec_mean_jrk = cat(2,mean_jrk(wts_range),mean_jrk(a2a_range));

[h_vel,sig_vel_mean]=ttest2(vec_mean_vel(:,1),vec_mean_vel(:,2));
[h_acc,sig_acc_mean]=ttest2(vec_mean_acc(:,1),vec_mean_acc(:,2));
[h_jrk,sig_jrk_mean]=ttest2(vec_mean_jrk(:,1),vec_mean_jrk(:,2));

%%
wd_bx = .8;
sz = 50;
figure()
subplot(131)
boxplot(vec_mean_vel,'Symbol', 'k.','Color','k','Widths',wd_bx);
hold on
for pos=1:size(vec_mean_vel,2)
    if ismember(pos,1:3:nxh1*3), c=1; else, c=2; end
    f=scatter(ones(1,max_n_group)+(pos-1)+(rand(01,max_n_group)-0.5)/10,vec_mean_vel(:,pos),sz,clrs_m(c,:),'filled','LineWidth',1.5);
end
hold off
%axis square;% ylim([-20 250])
XTickLabel={'CTR';'A2a'};
set(gca,'XTickLabel',XTickLabel); ylabel('velocity (px/s)');
set(gca,axeOpt{:})
title({'Mean velocity',sprintf('%s%1.3f','sig = ',sig_vel_mean)})


subplot(132)
boxplot(vec_mean_acc,'Symbol', 'k.','Color','k','Widths',wd_bx);
hold on
for pos=1:size(vec_mean_acc,2)
    if ismember(pos,1:3:nxh1*3), c=1; else, c=2; end
    f=scatter(ones(1,max_n_group)+(pos-1)+(rand(01,max_n_group)-0.5)/10,vec_mean_acc(:,pos),sz,clrs_m(c,:),'filled','LineWidth',1.5);
end
hold off
%axis square;% ylim([-20 250])
XTickLabel={'CTR';'A2a'};
set(gca,'XTickLabel',XTickLabel); ylabel('acceleration (px/s^2)');
set(gca,axeOpt{:})
title({'Mean acceleration',sprintf('%s%1.3f','sig = ',sig_acc_mean)})


subplot(133)
boxplot(vec_mean_jrk,'Symbol', 'k.','Color','k','Widths',wd_bx);
hold on
for pos=1:size(vec_mean_jrk,2)
    if ismember(pos,1:3:nxh1*3), c=1; else, c=2; end
    f=scatter(ones(1,max_n_group)+(pos-1)+(rand(01,max_n_group)-0.5)/10,vec_mean_jrk(:,pos),sz,clrs_m(c,:),'filled','LineWidth',1.5);
end
hold off
%axis square;% ylim([-20 250])
XTickLabel={'CTR';'A2a'};
set(gca,'XTickLabel',XTickLabel); ylabel('jerk (px/s^3)');
set(gca,axeOpt{:})
title({'Mean jerk',sprintf('%s%1.3f','sig = ',sig_jrk_mean)})

set(gcf,'position',[ 2323         222        1031         572])
saveas(gcf,strcat(save_path,filesep,'mean_kinematics.png'),'png')
saveas(gcf,strcat(save_path,filesep,'mean_kinematics.epsc'),'epsc')

%% Median vector 

med_vel = squeeze(median(median(speed_fwd,1,'omitnan'),2,'omitnan'));
med_acc = squeeze(median(median(acc_fwd,1,'omitnan'),2,'omitnan'));
med_jrk = squeeze(median(median(jerk_fwd,1,'omitnan'),2,'omitnan'));
nxh3 = 3;
vec_med_vel = cat(2,med_vel(wts_range),med_vel(a2a_range));
vec_med_acc = cat(2,med_acc(wts_range),med_acc(a2a_range));
vec_med_jrk = cat(2,med_jrk(wts_range),med_jrk(a2a_range));

[h_vel_med,sig_vel_med]=ttest2(vec_med_vel(:,1),vec_med_vel(:,2));
[h_acc_med,sig_acc_med]=ttest2(vec_med_acc(:,1),vec_med_acc(:,2));
[h_jrk_med,sig_jrk_med]=ttest2(vec_med_jrk(:,1),vec_med_jrk(:,2));

%%
wd_bx = .8;
sz = 50;
figure()
subplot(131)
boxplot(vec_med_vel,'Symbol', 'k.','Color','k','Widths',wd_bx);
hold on
for pos=1:size(vec_med_vel,2)
    if ismember(pos,1:3:nxh1*3), c=1; else, c=2; end
    f=scatter(ones(1,max_n_group)+(pos-1)+(rand(01,max_n_group)-0.5)/10,vec_med_vel(:,pos),sz,clrs_m(c,:),'filled','LineWidth',1.5);
end
hold off
%axis square;% ylim([-20 250])
XTickLabel={'CTR';'A2a'};
set(gca,'XTickLabel',XTickLabel); ylabel('velocity (px/s)');
set(gca,axeOpt{:})
title({'Median velocity',sprintf('%s%1.3f','sig = ',sig_vel_med)})


subplot(132)
boxplot(vec_med_acc,'Symbol', 'k.','Color','k','Widths',wd_bx);
hold on
for pos=1:size(vec_med_acc,2)
    if ismember(pos,1:3:nxh1*3), c=1; else, c=2; end
    f=scatter(ones(1,max_n_group)+(pos-1)+(rand(01,max_n_group)-0.5)/10,vec_med_acc(:,pos),sz,clrs_m(c,:),'filled','LineWidth',1.5);
end
hold off
%axis square;% ylim([-20 250])
XTickLabel={'CTR';'A2a'};
set(gca,'XTickLabel',XTickLabel); ylabel('acceleration (px/s^2)');
set(gca,axeOpt{:})
title({'Median acceleration',sprintf('%s%1.3f','sig = ',sig_acc_med)})


subplot(133)
boxplot(vec_med_jrk,'Symbol', 'k.','Color','k','Widths',wd_bx);
hold on
for pos=1:size(vec_med_jrk,2)
    if ismember(pos,1:3:nxh1*3), c=1; else, c=2; end
    f=scatter(ones(1,max_n_group)+(pos-1)+(rand(01,max_n_group)-0.5)/10,vec_med_jrk(:,pos),sz,clrs_m(c,:),'filled','LineWidth',1.5);
end
hold off
%axis square;% ylim([-20 250])
XTickLabel={'CTR';'A2a'};
set(gca,'XTickLabel',XTickLabel); ylabel('jerk (px/s^3)');
set(gca,axeOpt{:})
title({'Median jerk',sprintf('%s%1.3f','sig = ',sig_jrk_med)})

set(gcf,'position',[ 2323         222        1031         572])
saveas(gcf,strcat(save_path,filesep,'med_kinematics.png'),'png')
saveas(gcf,strcat(save_path,filesep,'med_kinematics.epsc'),'epsc')




