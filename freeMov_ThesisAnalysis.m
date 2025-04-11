%% Most important figures - fot the thesis
clear; close all; clc

%% % Paths for each mouse and initializaition variables
raw_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\raw_data';
mat_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\analyzed_data\mat_files';
out_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\analyzed_data\output_files';

%data_pathA = '/Users/teresaserradasduarte/Dropbox (Learning Lab)/Learning Lab Team Folder/Patlab protocols/data/TD/raw_data';
project_name = '20250106_A2aCasp_G2';
setup = 'freely_mov';
data_path = strcat(mat_folder,filesep,project_name,filesep,setup);
save_path = strcat(out_folder,filesep,project_name,filesep,setup,filesep,'group_results',filesep,'thesis');
if ~exist('save_path','dir'), mkdir(save_path); end
matname = 'group_results.mat';

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
fps = 120;
sess_dur_long = fps*60*45;

nr_trials_TOT = nan(n_sess_TOT,num_animals);
nr_trials_t = nan(n_sess,num_animals);
nr_reaches_t = nan(n_sess,num_animals);
length_sess = nan(n_sess,num_animals);
reach_times_sess = nan(nr_reaches_max_TOT,n_sess,num_animals);
dur_reaches_t = nan(nr_reaches_max,n_sess,num_animals);
dur_full_t = nan(nr_reaches_max,n_sess,num_animals);
r_paw_used_t = nan(nr_reaches_max,n_sess,num_animals);
start_fwd_t = nan(nr_reaches_max,n_sess,num_animals);
stop_reach_t = nan(nr_reaches_max,n_sess,num_animals);

cat_reach_id_t = nan(nr_reaches_max,n_sess,num_animals);
class_reaches_t = nan(3,nr_reaches_max,n_sess,num_animals);

iri_sec_t = nan(nr_reaches_max_TOT,n_sess,num_animals);
idx_iri_match_reach_t = nan(nr_reaches_max_TOT,n_sess,num_animals);
reach_xz_at_start_t = nan(2,nr_reaches_max_TOT,n_sess,num_animals);

time_to_1st_reach_peak_t = nan(nr_reaches_max,n_sess,num_animals);
time_to_1st_reach_start_t = nan(nr_reaches_max,n_sess,num_animals);
reaches_all_t = nan(nr_datapoints,2,nr_reaches_max,n_sess,num_animals);
lick_all_t = nan(nr_datapoints,2,nr_reaches_max,n_sess,num_animals);
reaches_start_stop_t = nan(nr_datapoints,2,nr_reaches_max,n_sess,num_animals);
speed_reaches_all_t = nan(nr_datapoints,2,nr_reaches_max,n_sess,num_animals);
distance_water_all_t = nan(nr_datapoints,nr_reaches_max,n_sess,num_animals);
speed_water_all_t = nan(nr_datapoints,nr_reaches_max,n_sess,num_animals);
water_pos_all_t = nan(2,n_sess,num_animals);

cat_reach_nr_t = nan(3,n_sess,num_animals);
puorpose_yesNo_nr = nan(2,n_sess,num_animals);
hit_yesNo_nr = nan(2,n_sess,num_animals);
suc_yesNo_nr = nan(2,n_sess,num_animals);

paw_ss_all = nan(sess_dur_long,2,n_sess,num_animals);


if ~exist(strcat(data_path,filesep,matname),"file")
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
            dur_full_t(1:nr_reaches_t(s,m),s,m) = reaches.start_stop.dur_start_stop_mat;
            r_paw_used_t(1:nr_reaches_t(s,m),s,m) = reaches.Rpaw_used;
            start_fwd_t(1:nr_reaches_t(s,m),s,m) = reaches.start_forw.start_forw_mat;
            stop_reach_t(1:nr_reaches_t(s,m),s,m) = reaches.start_stop.stop_reach_mat;
            paw_ss_all(1:size(sessions.paw_ss,1),:,s,m) = sessions.paw_ss;
            lick_all_t(:,:,1:nr_reaches_t(s,m),s,m) = reaches.lick_reach_mat;


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
                'reach_to_water','speed_reach_to_water','water_pos',...
                'perc1','perc2','perc3','porpuse_nr','noporpuse_nr',...
                'hit_nr','miss_nr','succ_nr','fail_nr',...
                'idx_match_iri','reach_xz_at_start');

            iri_sec_t(1:length(iri_sec),s,m) = iri_sec;
            idx_iri_match_reach_t(1:length(idx_match_iri),s,m) = idx_match_iri;
            reach_xz_at_start_t(:,1:size(reach_xz_at_start,2),s,m) = reach_xz_at_start;

            time_to_1st_reach_peak_t(1:length(time_to_1st_reach_peak),s,m) = time_to_1st_reach_peak;
            time_to_1st_reach_start_t(1:length(time_to_1st_reach_start),s,m) = time_to_1st_reach_start;
            reaches_all_t(:,:,1:nr_reaches_t(s,m),s,m) = reaches_all;
            reaches_start_stop_t(:,:,1:nr_reaches_t(s,m),s,m) = reaches_start_stop;
            speed_reaches_all_t(:,:,1:nr_reaches_t(s,m),s,m) = speed_reaches_all;
            distance_water_all_t(:,1:nr_reaches_t(s,m),s,m) = reach_to_water;
            speed_water_all_t(:,1:nr_reaches_t(s,m),s,m) = speed_reach_to_water;
            water_pos_all_t(:,s,m) = water_pos';

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
fprintf(['%s' '\n'],'saving...')
save(strcat(data_path,filesep,matname))
fprintf(['%s' '\n'],'done!')

else
    load(strcat(data_path,filesep,matname))
    fprintf(['%s' ...
    '\n'],'loaded!')
end


%% Convert reaches from px to mm
fm_px2mm = .1;
reaches_mm = reaches_all_t.*fm_px2mm;
reaches_mm_allSess = paw_ss_all.*fm_px2mm;
licks_mm = lick_all_t.*fm_px2mm;

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

lw = 1.5;
sessions_TOT = {'S1';'S2';'S3';'S4';'S5'};
axeOpt = {'linewidth',1.5,'box','off','GridAlpha',0.05,'ticklength',[1,1]*.01,'TickDir','out'};
max_n_group = max(length(wts_range),length(a2a_range));
nxh1 = 1;
wd_bx = .8;
sz = 50;

%% LEARNING - for thesis
disp_nsess = 3;
trans = .15;
lw_ind = 1;
lw_m = 2;

subplot(221)
plot(nr_reaches_t(:,wts_range),'-','LineWidth',lw_ind,'Color',cat(2,wts_clr,trans)); hold on
plot(nr_reaches_t(:,a2a_range),'-','LineWidth',lw_ind,'Color',cat(2,a2a_clr,trans));
a=plot(mean(nr_reaches_t(:,wts_range),2),'-','LineWidth',lw_m,'Color',wts_clr); hold on
b=plot(mean(nr_reaches_t(:,a2a_range),2),'-','LineWidth',lw_m,'Color',a2a_clr); hold on
ylim([0 1300])
hold off

set(gca,axeOpt{:})
%title('Nr of reaches per session')
xlim([.9 n_sess+.1])
%ylim([280 1100])
xticks(1:n_sess)
xticklabels(sessions_name)
xlabel('sessions'); ylabel('number of reaches')
legend([a,b],'CTR','A2a','box','off','location','southeast')
%legend(cat(1,mice(wts_range),mice(a2a_range)),'box','off','location','southeastoutside','Interpreter', 'none')
%saveas(gcf,strcat(save_path,filesep,'number_reaches.png'),'png')

subplot(222)
plot(nr_trials_TOT(:,wts_range),'-','LineWidth',lw_ind,'Color',cat(2,wts_clr,trans)); hold on
plot(nr_trials_TOT(:,a2a_range),'-','LineWidth',lw_ind,'Color',cat(2,a2a_clr,trans));
plot(mean(nr_trials_TOT(:,wts_range),2),'-','LineWidth',lw_m,'Color',wts_clr); 
plot(mean(nr_trials_TOT(:,a2a_range),2),'-','LineWidth',lw_m,'Color',a2a_clr); hold off
xlim([.8 n_sess_TOT+.2])
ylim([50 400])
set(gca,axeOpt{:})
%title('Nr of trials per session')
xticks(1:n_sess_TOT)
xticklabels(sessions_TOT)
xlabel('sessions'); ylabel('number of trials')
%saveas(gcf,strcat(save_path,filesep,'number_trials_session_learning.png'),'png')


subplot(223)
nr_suc = squeeze(suc_yesNo_nr(1,:,:));

plot(nr_suc(:,wts_range),'-','LineWidth',lw_ind,'Color',cat(2,wts_clr,trans)); hold on
plot(nr_suc(:,a2a_range),'-','LineWidth',lw_ind,'Color',cat(2,a2a_clr,trans));
plot(mean(nr_suc(:,wts_range),2),'-','LineWidth',lw_m,'Color',wts_clr); 
plot(mean(nr_suc(:,a2a_range),2),'-','LineWidth',lw_m,'Color',a2a_clr); hold off

set(gca,axeOpt{:})
%title('success reaches')
xlim([disp_nsess-.3 n_sess+.1])
xticks(disp_nsess:n_sess)
xticklabels(sessions_name(disp_nsess:n_sess))
xlabel('sessions'); ylabel('success reaches')
%legend(cat(1,mice(wts_range),mice(a2a_range)),'box','off','location','southeastoutside','Interpreter', 'none')

 last_suc = mean(nr_suc(4:5,:),1);
 vec_nr_suc =nan(max_n_group,(nxh1*3)-1);
 vec_nr_suc(1:length(wts_range),1:3:nxh1*3)=last_suc(wts_range)';
 vec_nr_suc(1:length(a2a_range),2:3:nxh1*3+1)=last_suc(a2a_range)';
 length_n_suc=size(vec_nr_suc,1);
 [h_nsuc,sig_nsuc]=ttest2(vec_nr_suc(:,1),vec_nr_suc(:,2));

subplot(224)
boxplot(vec_nr_suc,'Symbol', 'k.','Color','k','Widths',0.8);
hold on
for pos=1:size(vec_nr_suc,2)
    if ismember(pos,1:3:nxh1*3), c=1; else, c=2; end
    f=scatter(ones(1,length_n_suc)+(pos-1).*(1+(rand(01,length_n_suc)-0.5)/50),vec_nr_suc(:,pos),[],clrs_m(c,:),'filled','LineWidth',1.5);
end
hold off
axis square;% ylim([-20 250])
XTickLabel={'CTR';'DTX'};
set(gca,'XTickLabel',XTickLabel); ylabel('Nr of trials');
set(gca,'linewidth',1.2,'box','off','layer','top','GridAlpha',0.05);
title({'Nr successes',sprintf('%s%1.3f','sig = ',sig_nsuc)})

%saveas(gcf,strcat(save_path,filesep,'percent_suc_over_purpose.png'),'png')
set(groot, 'DefaultAxesFontName', 'Arial', 'DefaultTextFontName', 'Arial');
set(gcf,'Position',[2453         110         980         832], 'Color','w')
saveas(gcf,strcat(save_path,filesep,'learning_thesis.png'),'png')
saveas(gcf,strcat(save_path,filesep,'learning_thesis.epsc'),'epsc')


%% Speed forward only
reaches_forward_part = nan(41,2,nr_reaches_max,n_sess,num_animals);

for m=1:num_animals
    for s=1:n_sess
        for rr = 1:nr_reaches_max
            if ~isnan(start_fwd_t(rr,s,m))
                reaches_forward_part(start_fwd_t(rr,s,m):41,:,rr,s,m) = reaches_mm(start_fwd_t(rr,s,m):41,:,rr,s,m);
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

% Instantaneous acelleration
acc_fwd_partial = nan(size_npt-2,2,nr_reaches_max,num_animals);
for i=1:size_npt-2
    acc_fwd_partial(i,:,:,:) = reaches_fwd_last(i+2,:,:,:) - 2.*reaches_fwd_last(i+1,:,:,:) + reaches_fwd_last(i,:,:,:);
end
acc_fwd = squeeze(sqrt(acc_fwd_partial(:,1,:,:).^2+acc_fwd_partial(:,2,:,:).^2)./dt.^2);

% Instantaneous jerk
jerk_fwd_partial = nan(size_npt-3,2,nr_reaches_max,num_animals);
for i=1:size_npt-3
    jerk_fwd_partial(i,:,:,:) = reaches_fwd_last(i+3,:,:,:) - 3*reaches_fwd_last(i+2,:,:,:) ...
        + 3.*reaches_fwd_last(i+1,:,:,:) - reaches_fwd_last(i,:,:,:);
end
jerk_fwd = squeeze(sqrt(jerk_fwd_partial(:,1,:,:).^2+jerk_fwd_partial(:,2,:,:).^2)./dt.^3);


% Median vector 
med_vel = squeeze(median(median(speed_fwd,1,'omitnan'),2,'omitnan'));
med_acc = squeeze(median(median(acc_fwd,1,'omitnan'),2,'omitnan'));
med_jrk = squeeze(median(median(jerk_fwd,1,'omitnan'),2,'omitnan'));
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
set(gca,'XTickLabel',XTickLabel); ylabel('velocity (mm/s)');
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
set(gca,'XTickLabel',XTickLabel); ylabel('acceleration (mm/s^2)');
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
set(gca,'XTickLabel',XTickLabel); ylabel('jerk (mm/s^3)');
set(gca,axeOpt{:})
title({'Median jerk',sprintf('%s%1.3f','sig = ',sig_jrk_med)})

set(gcf,'position',[ 2323         222        1031         572])
saveas(gcf,strcat(save_path,filesep,'med_kinematics.png'),'png')
saveas(gcf,strcat(save_path,filesep,'med_kinematics.epsc'),'epsc')


%%
% Tranjectories examples
selected_represent = [1, 7, 4, 2, 9, 10];

%reaches_rr_all(reaches_rr_all==0)=nan;
transpa = 0.03;
x_lim = [0 350].*fm_px2mm;
z_lim = [0 450].*fm_px2mm;
s=5;
%for s=1:num_sessions
    figure();
    tt = tiledlayout(2,3);
    title(tt,sprintf('%s%s','reach trajectories ',char(sessions_name(s))),'Interpreter','none')
    for m = selected_represent
        nexttile
        plot(squeeze(reaches_mm(:,1,:,s,m)),...
            squeeze(reaches_mm(:,2,:,s,m)),"Color",cat(2,clrs_m(m,:),transpa));
        hold on
%                  plot(mean(reaches_all_t(:,1,:,s,m),3,'omitnan'),...
%                      mean(reaches_all_t(:,2,:,s,m),3,'omitnan'),...
%                      "Color",clrs_m(m,:),'LineWidth',2);
        set(gca,axeOpt{:})
        axis square;
        xlabel('x'); ylabel('y')
        axis([x_lim z_lim])
        title(char(mice(m)),'Interpreter','none')
    end
    set(gcf,'Color','w', 'position',[1975         161        1203         821])

saveas(gcf,strcat(save_path,filesep,'traject_representatives.png'),'png')
saveas(gcf,strcat(save_path,filesep,'traject_representatives.epsc'),'epsc')
print(gcf, fullfile(save_path, 'traject_representatives.pdf'), '-dpdf', '-painters');


%% Speed towards water

water_all_reps = permute(repmat(water_pos_all_t(:,:,:).*fm_px2mm,[1,1,1,size(reaches_mm,1),size(reaches_mm,3)]),[4 1 5 2 3]);
distance_to_water_mm = squeeze(sqrt(...
            (reaches_mm(:,1,:,:,:)-water_all_reps(:,1,:,:,:)).^2 + ...
            (reaches_mm(:,2,:,:,:)-water_all_reps(:,2,:,:,:)).^2 ...
            ));
speed_to_water_mm = cat(1,nan(1,nr_reaches_max,n_sess,num_animals),diff(distance_to_water_mm,1));



%%
transpaa = .5;
sess=5;
lw_m = 1.5;
lw_mm = 3;
figure()
% subplot(211)
% plot(tm,squeeze(median(distance_to_water_mm(:,:,sess,wts_range),2,'omitnan')),'linewidth',lw_m,'color',cat(2,wts_clr,transpaa)); hold on
% plot(tm,squeeze(median(distance_to_water_mm(:,:,sess,a2a_range),2,'omitnan')),'linewidth',lw_m,'color',cat(2,a2a_clr,transpaa)); hold on
% 
% plot(tm,mean(squeeze(median(distance_to_water_mm(:,:,sess,wts_range),2,'omitnan')),2),'linewidth',lw_mm,'color',wts_clr); hold on
% plot(tm,mean(squeeze(median(distance_to_water_mm(:,:,sess,a2a_range),2,'omitnan')),2),'linewidth',lw_mm,'color',a2a_clr); hold on
% xlabel('time (s)'); ylabel('speed towards water (px/s)')
% xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
% set(gca,axeOpt{:})
% 
% subplot(212)
plot(tm,squeeze(median(speed_to_water_mm(:,:,sess,wts_range),2,'omitnan')),'linewidth',lw_m,'color',cat(2,wts_clr,transpaa)); hold on
plot(tm,squeeze(median(speed_to_water_mm(:,:,sess,a2a_range),2,'omitnan')),'linewidth',lw_m,'color',cat(2,a2a_clr,transpaa)); hold on

a=plot(tm,mean(squeeze(median(speed_to_water_mm(:,:,sess,wts_range),2,'omitnan')),2),'linewidth',lw_mm,'color',wts_clr); hold on
b=plot(tm,mean(squeeze(median(speed_to_water_mm(:,:,sess,a2a_range),2,'omitnan')),2),'linewidth',lw_mm,'color',a2a_clr); hold on
xlabel('time (s)'); ylabel('speed towards water (mm/s)')
xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
set(gca,axeOpt{:})

set(gcf,'position',[2488         362         968         407],'Color','w');
legend([a,b],'CTR','A2a','box','off')
saveas(gcf,strcat(save_path,filesep,'speed2water.png'),'png')
saveas(gcf,strcat(save_path,filesep,'speed2water.epsc'),'epsc')
print(gcf, fullfile(save_path, 'speed2water.pdf'), '-dpdf', '-painters');

% %%
% transpa=0.03;
% m=4
% figure()
% plot(tm,speed_water_all_t(:,:,sess,m),'linewidth',lw_m,'color',cat(2,clrs_m(m,:),transpa)); hold on
% hold on
% plot(tm,mean(squeeze(median(speed_to_water_mm(:,:,sess,m),2,'omitnan')),2),'linewidth',lw_mm,'color',clrs_m(m,:)); hold on
% xlabel('time (s)'); ylabel('speed towards water (px/s)')
% xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
% set(gcf,'position',[2488         362         968         407],'Color','w');
% set(gca,axeOpt{:})

%% Distance to mean - distributions
% SPACE-ONLY -------------------------
   % Interplolate each single trial
   max_reach = 41;
   timepoints = 1:max_reach;
   n_timepoints = length(timepoints);
   n_dims =  size(reaches_mm,2);
   res_or = timepoints;
   res_high = 0.1:0.1:n_timepoints;
   reaches_mm_mean = squeeze(mean(reaches_mm(timepoints,:,:,:,:),3,'omitnan'));
   reaches_mm_med = squeeze(median(reaches_mm(timepoints,:,:,:,:),3,'omitnan'));
   nbinsplus=21;
    nbins = nbinsplus-1;

   % Allocatre space
   % Highres trial trajectories
   hires_reach_m = nan(numel(res_high),n_dims,nr_reaches_max,n_sess,num_animals);
   % distance
   reach_distMed_sp_m = nan(n_timepoints,nr_reaches_max,n_sess,num_animals);
   % Integral of distance
   reach_distMed_sp_int_m = nan(nr_reaches_max,n_sess,num_animals);
   % bins and cdf
   binsMed_reach_sp_m = nan(nbinsplus,n_sess,num_animals);
   binsMed_reach_sp_norm_m = nan(nbinsplus,n_sess,num_animals);
   cdfMean_reach_sp_norm_m = nan(nbinsplus,n_sess,num_animals);

   for m =1:num_animals
       for s = 1:n_sess
           for t = 1:nr_reaches_max
               if sum(isnan(reaches_mm(timepoints,1,t,s,m)))<n_timepoints
                   for dim = 1:n_dims
                       tmp_xy = reaches_mm(timepoints,dim,t,s,m);
                       last_nan =find(isnan(tmp_xy),1,'last')+1;
                       tmp_xy(isnan(tmp_xy))=0;
                       hi_tmp_xy = spline(res_or,tmp_xy,res_high);
                       hi_tmp_xy(1:last_nan*10)=nan;
                       hires_reach_m(:,dim,t,s,m) = hi_tmp_xy;
                   end
                   for pt=1:n_timepoints
                       % DISTANCE TO  THE MED
                       tmp_dist_allMed=sqrt((hires_reach_m(:,1,t,s,m)-reaches_mm_med(pt,1,s,m)).^2+...
                           (hires_reach_m(:,2,t,s,m)-reaches_mm_med(pt,2,s,m)).^2);
                       [min_dist_med,min_dist_loc_med]=min(tmp_dist_allMed);
                       reach_distMed_sp_m(pt,t,s,m)=min_dist_med;
                   end

               end
           end

           %%

              range_reach = 30:41;
           % distribution
    reach_distMed_sp_int_m(:,s,m) = squeeze(sum(reach_distMed_sp_m(range_reach,:,s,m),1)); % integral
    % cumulative distribution function
    [countsMed_reach_all_sp, binsMed_reach_sp_m(:,s,m)] = histcounts(reach_distMed_sp_int_m(:,s,m),nbins);
    cdfMed_reach_all_sp = cumsum(countsMed_reach_all_sp); cdfMed_reach_all_sp(end+1)=cdfMed_reach_all_sp(end);
    cdfMean_reach_sp_norm_m(:,s,m)=cdfMed_reach_all_sp./cdfMed_reach_all_sp(end);
       end
   end

   % Last session pool
   % ctr
   [counts_ctr, bins_distMed_ctr] = histcounts(reach_distMed_sp_int_m(:,5,wts_range),nbins);
   cdfMed_ctr = cumsum(counts_ctr); cdfMed_ctr(end+1)=cdfMed_ctr(end);
   cdfMed_ctr_norm=cdfMed_ctr./cdfMed_ctr(end);
    %a2a
   [counts_a2a, bins_distMed_a2a] = histcounts(reach_distMed_sp_int_m(:,5,a2a_range),nbins);
   cdfMed_a2a = cumsum(counts_a2a); cdfMed_a2a(end+1)=cdfMed_a2a(end);
   cdfMed_a2a_norm=cdfMed_a2a./cdfMed_a2a(end);

  %%
m=5
figure(1)
subplot(121)
plot(reaches_mm_med(range_reach,1,5,m),reaches_mm_med(range_reach,2,5,m),'linewidth',lw_m,'color',cat(2,clrs_m(m,:),1))
subplot(122)
plot(reach_distMed_sp_m(range_reach,:,5,m),'linewidth',lw_m,'color',cat(2,clrs_m(m,:),transpa)); hold on
plot(median(reach_distMed_sp_m(range_reach,:,5,m),2,'omitnan'),'linewidth',lw_m,'color',cat(2,clrs_m(m,:),1))


   %% Vector of med
% Median vector 
med_distMed = squeeze(median(reach_distMed_sp_int_m(:,5,:),1,'omitnan'));
vec_med_distMed = cat(2,med_distMed(wts_range),med_distMed(a2a_range));
[h_distMed,sig_distMed] = ttest2(vec_med_distMed(:,1),vec_med_distMed(:,2));


   %%
   figure
   for m = 1:num_animals
       stairs(binsMed_reach_sp_m(:,5,m),cdfMean_reach_sp_norm_m(:,5,m),'color',cat(2,clrs_m(m,:),transpa),'linewidth', .5); hold on
   end
   stairs(bins_distMed_ctr,cdfMed_ctr_norm,'color',cat(2,wts_clr,1),'linewidth', 2); hold on
   stairs(bins_distMed_a2a,cdfMed_a2a_norm,'color',cat(2,a2a_clr,1),'linewidth', 2); hold off
   xlim([0 100])
   xlabel('sum distance to median'); ylabel('cdf')
   set(gca,axeOpt{:})
saveas(gcf,strcat(save_path,filesep,'cdf_variance.png'),'png')
saveas(gcf,strcat(save_path,filesep,'cdf_variance.epsc'),'epsc')

   %
   figure
  wd_bx = .8;
boxplot(vec_med_distMed,'Symbol', 'k.','Color','k','Widths',wd_bx);
hold on
for pos=1:size(vec_med_distMed,2)
    if ismember(pos,1:3:nxh1*3), c=1; else, c=2; end
    f=scatter(ones(1,max_n_group)+(pos-1)+(rand(01,max_n_group)-0.5)/10,vec_med_distMed(:,pos),sz,clrs_m(c,:),'filled','LineWidth',1.5);
end
hold off
%axis square;% ylim([-20 250])
XTickLabel={'CTR';'A2a'};

set(gca,'XTickLabel',XTickLabel); ylabel('distance to median');
set(gca,axeOpt{:})
title({'Median variance',sprintf('%s%1.3f','sig = ',sig_distMed)})
set(gcf,'position',[ 3341         239         179         391]);
saveas(gcf,strcat(save_path,filesep,'cdf_variance_inset.png'),'png')
saveas(gcf,strcat(save_path,filesep,'cdf_variance_inset.epsc'),'epsc')

%% END-POINT IN X AND Z
mean_pos_spout=[227,240]*fm_px2mm;
figure
figg=tiledlayout(1,2);
bw = .4;
%axis_xx = [190 350];
range_histogram = 1:61;
sesses = 5;
sess_interval_name = sprintf('%s%i%s%i','S',sesses(1),'-',sesses(end));

nexttile
end_point_x_wt = max(reaches_mm(range_histogram,1,:,sesses,wts_range));
end_point_x_a2 = max(reaches_mm(range_histogram,1,:,sesses,a2a_range));
h1=histogram(end_point_x_wt(~isnan(end_point_x_wt)),'binwidth',bw,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
h2=histogram(end_point_x_a2(~isnan(end_point_x_a2)),'binwidth',bw,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold on
xline(mean_pos_spout(1),'--','LineWidth',2,'color',[.9 .9 .9])
%xlim(axis_xx)
set(gca,axeOpt{:})
xlabel('end point in x (mm)'); ylabel('pdf')
legend('wts','a2as','box','off')

nexttile
bw = .8;
end_point_z_wt = max(reaches_mm(range_histogram,2,:,sesses,wts_range));
end_point_z_a2 = max(reaches_mm(range_histogram,2,:,sesses,a2a_range));
h1=histogram(end_point_z_wt(~isnan(end_point_z_wt)),'binwidth',bw,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
h2=histogram(end_point_z_a2(~isnan(end_point_z_a2)),'binwidth',bw,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold on
xline(mean_pos_spout(2),'--','LineWidth',2,'color',[.9 .9 .9])
%xlim([150 450])
set(gca,axeOpt{:})
xlabel('end point in z (mm)'); ylabel('pdf')
legend('wts','a2as','box','off')

set(gcf,'color','w','position',[2646         351         874         368])
title(figg,sprintf('%s%s%s%s','Distribution of end points pooled ' ,...
    char(sessions_name(sesses(1))),'-',char(sessions_name(sesses(end)))))

% saveas(gcf,strcat(save_path,filesep,'endpoints_pooled_',...
%     char(sessions_name(sesses(1))),'-',char(sessions_name(sesses(end)))),'png')
% 
print(gcf, fullfile(save_path, 'endpoint.pdf'), '-dpdf', '-painters');

%%
% Median vector 
range_histogram = 1:61;
nxh1 = 1;
mean_endpointX = squeeze(mean(max(reaches_mm(range_histogram,1,:,sesses,:)),3,'omitnan'));
mean_endpointZ =  squeeze(mean(max(reaches_mm(range_histogram,2,:,sesses,:)),3,'omitnan'));
vec_mean_epX = cat(2,mean_endpointX(wts_range),mean_endpointX(a2a_range));
vec_mean_epZ = cat(2,mean_endpointZ(wts_range),mean_endpointZ(a2a_range));

[h_epX_mean,sig_epX_mean]=ttest2(vec_mean_epX(:,1),vec_mean_epX(:,2));
[h_epZ_mean,sig_epZ_mean]=ttest2(vec_mean_epZ(:,1),vec_mean_epZ(:,2));


  wd_bx = .4;

   figure
   subplot(121)
boxplot(vec_mean_epX,'Symbol', 'k.','Color','k','Widths',wd_bx);
hold on
for pos=1:size(vec_mean_epX,2)
    if ismember(pos,1:3:nxh1*3), c=1; else, c=2; end
    f=scatter(ones(1,max_n_group)+(pos-1)+(rand(01,max_n_group)-0.5)/10,vec_mean_epX(:,pos),sz,clrs_m(c,:),'filled','LineWidth',1.5);
end
hold off
%axis square;% ylim([-20 250])
XTickLabel={'CTR';'A2a'};

set(gca,'XTickLabel',XTickLabel); ylabel('mean end-point in x (mm)');
set(gca,axeOpt{:})
title({'Mean end-point in x',sprintf('%s%1.3f','sig = ',sig_epX_mean)})
%set(gcf,'position',[ 3341         239         179         391]);

   subplot(122)
boxplot(vec_mean_epZ,'Symbol', 'k.','Color','k','Widths',wd_bx);
hold on
for pos=1:size(vec_mean_epZ,2)
    if ismember(pos,1:3:nxh1*3), c=1; else, c=2; end
    f=scatter(ones(1,max_n_group)+(pos-1)+(rand(01,max_n_group)-0.5)/10,vec_mean_epZ(:,pos),sz,clrs_m(c,:),'filled','LineWidth',1.5);
end
hold off
%axis square;% ylim([-20 250])
XTickLabel={'CTR';'A2a'};

set(gca,'XTickLabel',XTickLabel); ylabel('mean end-point in z (mm)');
set(gca,axeOpt{:})
title({'Mean end-point in z',sprintf('%s%1.3f','sig = ',sig_epZ_mean)})
%set(gcf,'position',[ 3341         239         179         391]);

% saveas(gcf,strcat(save_path,filesep,'cdf_variance_inset.png'),'png')
% saveas(gcf,strcat(save_path,filesep,'cdf_variance_inset.epsc'),'epsc')
print(gcf, fullfile(save_path, 'inset_endpoint.pdf'), '-dpdf', '-painters');




% ----------------------------------------
%% iri vs Position of start - individual animals
figure
bw=.2;
axis_xx = [0 2];
for sess = 1:num_sessions
    figg=tiledlayout(2,num_animals/2);
   
    for m = 1:num_animals
        nexttile
        
        scatter(iri_sec_t(:,sess,m),reach_xz_at_start_t(1,:,sess,m), ...
            10,clrs_m(m,:),'filled'); hold on
        xlim(axis_xx)
        %ylim([50 300])
        set(gca,axeOpt{:})
        xlabel('iri (s)'); ylabel('position (px)')
        title(mice(m),'color',clrs_m(m,:),'Interpreter','none');
    end
    title(figg,{'inter-reach interval';char(sessions_name(sess))})
    set(gcf,'position',[1950         159        1830         729],'color','w')
    %saveas(gcf,strcat(save_path,filesep,'iri',char(sessions_name(sess)),'.png'),'png')

end

%%
sess = 5;
iri_pool_a2 = squeeze(iri_sec_t(2:end,sess,a2a_range));
iri_pool_wt = squeeze(iri_sec_t(2:end,sess,wts_range));
reachX_at_start_wts = squeeze(reach_xz_at_start_t(1,2:end,sess,wts_range));
reachX_at_start_a2a = squeeze(reach_xz_at_start_t(1,2:end,sess,a2a_range)); 
mask_nonan_wts = ~isnan(reachX_at_start_wts);
mask_nonan_a2a = ~isnan(reachX_at_start_a2a);

iri_wts = iri_pool_wt(mask_nonan_wts);
iri_a2a = iri_pool_a2(mask_nonan_a2a);
startX_wt = reachX_at_start_wts(mask_nonan_wts).*fm_px2mm;
startX_a2 = reachX_at_start_a2a(mask_nonan_a2a).*fm_px2mm;


%% iri
figure
figg=tiledlayout(1,1);
sesses = 5;
bw=.1;
axis_xx = [-.5 6];
iri_lims = .85;


h1=histogram(iri_wts,'binwidth',bw,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
h2=histogram(iri_a2a,'binwidth',bw,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold on
xlim(axis_xx)
xline(iri_lims,'--','LineWidth',2,'color',[.9 .9 .9])
xline(dur_full_wts_med,'--','LineWidth',1.5,'color',wts_clr,'Alpha',0.8)
xline(dur_full_a2a_med,'--','LineWidth',1.5,'color',a2a_clr,'Alpha',0.5)
xline(quantile(iri_wts,.25),'-','LineWidth',1,'color',wts_clr,'Alpha',0.2)
xline(quantile(iri_a2a,.25),'-','LineWidth',1,'color',a2a_clr,'Alpha',0.2)

set(gca,axeOpt{:})

xlabel('iri (s)'); ylabel('pdf')
set(gcf,'color','w','position',[2493         275         603         499])
sess_interval_name = sprintf('%s%i%s%i','S',sesses(1),'-',sesses(end));
title(figg,sprintf('%s%i%s%i','Inter-reach interval S',sesses(1),'-',sesses(end)))

%saveas(gcf,strcat(save_path,filesep,'iri_pooled_,',sess_interval_name,'.png'),'png')
print(gcf, fullfile(save_path, 'iri_pool.pdf'), '-dpdf', '-painters');

%% MED IRI
  wd_bx = .4;

  % Median vector 
nxh1 = 1;
med_iri = squeeze(median(iri_sec_t(:,sesses,:),1,'omitnan'));
vec_med_iri = cat(2,med_iri(wts_range),med_iri(a2a_range));
[h_iri,sig_iri]=ttest2(vec_med_iri(:,1),vec_med_iri(:,2));


   figure
boxplot(vec_med_iri,'Symbol', 'k.','Color','k','Widths',wd_bx);
hold on
for pos=1:size(vec_med_iri,2)
    if ismember(pos,1:3:nxh1*3), c=1; else, c=2; end
    f=scatter(ones(1,max_n_group)+(pos-1)+(rand(01,max_n_group)-0.5)/10,vec_med_iri(:,pos),sz,clrs_m(c,:),'filled','LineWidth',1.5);
end
hold off
%axis square;% ylim([-20 250])
XTickLabel={'CTR';'A2a'};

set(gca,'XTickLabel',XTickLabel); ylabel('median IRI (s)');
set(gca,axeOpt{:})
title({'Mean end-point in x',sprintf('%s%1.3f','sig = ',sig_iri)})
set(gcf,'position',[ 2653         452         182         420]);
print(gcf, fullfile(save_path, 'iri_med_vec.pdf'), '-dpdf', '-painters');

%% TOGETHER
lim_iri = [0 6];
lim_posX = [2 20];

figure
% Create a tiled layout (2 rows, 2 cols) with tight spacing
t = tiledlayout(3,3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Top histogram (aligned with scatter plot)
ax1 = nexttile(6,[2,1]);
bw=.5;
histogram(startX_wt,'binwidth',bw,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5, 'Orientation', 'horizontal'); hold on
histogram(startX_a2,'binwidth',bw,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5, 'Orientation', 'horizontal'); hold on
ylim(lim_posX); 
set(ax1, 'XTick', [], 'YTick', []); % Remove ticks
ax1.XAxis.Visible = 'off'; % Hide X-axis
ax1.YAxis.Visible = 'off'; % Hide X-axis
%set(ax1,'XDir','reverse')
%xlim(lim_posX); 
%ylabel('pdf')

% Right-side histogram (aligned with scatter plot)
ax2 = nexttile(1,[1,2]);
 bw=.05;
histogram(iri_wts,'binwidth',bw,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
histogram(iri_a2a,'binwidth',bw,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold on
set(ax2, 'XTick', [], 'YTick', []); % Remove ticks
%set(ax2,'YDir','reverse')
ax2.YAxis.Visible = 'off'; % Hide Y-axis
ax2.XAxis.Visible = 'off'; % Hide X-axis
xlim(lim_iri); 
 set(gca, 'XScale', 'log');

szz = 15;
% Central scatter plot (spanning 2x2 grid)
ax3 = nexttile(4, [2,2]); % Scatter plot in center
 scatter(iri_wts,startX_wt, ...
            szz,wts_clr,'filled','MarkerFaceAlpha',.8);
hold on
  scatter(iri_a2a,startX_a2, ...
            szz,a2a_clr,'filled','MarkerFaceAlpha',.5);
 xlim(lim_iri)
 set(gca, 'XScale', 'log');
ylim(lim_posX)
 set(gca,axeOpt{:})
 xlabel('iri (s)'); ylabel('x at reach start (mm)')

set(gcf,'position',[2515 312 701 614],'color','w')

print(gcf, fullfile(save_path, 'scatter_startPosX_IRI.pdf'), '-dpdf', '-painters');

 %% duration and quantile

 figure
subplot(211)
bw=.02;
 histogram(dur_full_wts(~isnan(dur_full_wts)),'binwidth',bw,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
histogram(dur_full_a2a(~isnan(dur_full_a2a)),'binwidth',bw,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold on
%xlim(axis_xx)
%xline(iri_lims,'--','LineWidth',2,'color',[.9 .9 .9])
xline(dur_full_wts_med,'--','LineWidth',2,'color',wts_clr,'Alpha',0.8)
xline(dur_full_a2a_med,'--','LineWidth',2,'color',a2a_clr,'Alpha',0.5)
xline(iri_lims,'--','LineWidth',2,'color',[.9 .9 .9])
set(gca,axeOpt{:})
xlabel('Duration of full reach (s)'); ylabel('pdf')

 subplot(212)
idx_q25_iri_a2a = find(iri_a2a<quantile(iri_a2a,.25));
idx_q25_iri_wts = find(iri_wts<quantile(iri_wts,.25));
bw=.5;
histogram(startX_wt(idx_q25_iri_wts),'binwidth',bw,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
histogram(startX_a2(idx_q25_iri_a2a),'binwidth',bw,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold on
set(gca,axeOpt{:})
xlabel('x at reach start for short iri (mm)'); ylabel('pdf')
xlim([2 20])
set(gcf,'color','w','position',[2493         275         603         499])

 print(gcf, fullfile(save_path, 'duration_and_quantile.pdf'), '-dpdf', '-painters');


 
 %% Check all session

time_sess = (0:size(reaches_mm_allSess,1)-1)./fps;
time_selec =  3.5*60*fps:4*60*fps;
x_lim = [80 300].*fm_px2mm;
z_lim = [100 370].*fm_px2mm;
pos_rec_x = [time_sess(time_selec(1))./60 x_lim(1) ...
    (time_sess(time_selec(end))./60)-(time_sess(time_selec(1))./60)...
    x_lim(end)-x_lim(1)];
pos_rec_z = [time_sess(time_selec(1))./60 z_lim(1) ...
    (time_sess(time_selec(end))./60)-(time_sess(time_selec(1))./60)...
    z_lim(end)-z_lim(1)];

selected_represent = [7,2];

s=5;
for m=selected_represent
    figure();
    tt = tiledlayout(4,1);

    title(tt,sprintf('%s%s%s',char(mice(m)), ' , ',char(sessions_name(s))),'Interpreter','none')

    nexttile
    plot(time_sess./60,reaches_mm_allSess(:,1,s,m),'color',clrs_m(m,:));
    rectangle('Position',pos_rec_x, 'EdgeColor',[.8 .8 .8],'LineWidth',2)
    set(gca,axeOpt{:})
    ylim([x_lim])
    xlim([0 25])
    xlabel('time in session (min)');
    ylabel('x (px)');

    nexttile
    plot(time_sess./60,reaches_mm_allSess(:,2,s,m),'color',clrs_m(m,:));
    rectangle('Position',pos_rec_z, 'EdgeColor',[.8 .8 .8],'LineWidth',2)
    set(gca,axeOpt{:})
    ylim([z_lim])
    xlim([0 30])
    xlabel('time in session (min)');
    ylabel('z (px)');


    nexttile
    plot(time_sess(time_selec)./60,reaches_mm_allSess(time_selec,1,s,m),'color',clrs_m(m,:));
    set(gca,axeOpt{:})
    ylim([x_lim])
    xlim([time_sess(time_selec(1))./60 time_sess(time_selec(end))./60])
    xlabel('time (min)');
    ylabel('x (px)');

    nexttile
    plot(time_sess(time_selec)./60,reaches_mm_allSess(time_selec,2,s,m),'color',clrs_m(m,:));
    set(gca,axeOpt{:})
    ylim([z_lim])
    xlim([time_sess(time_selec(1))./60 time_sess(time_selec(end))./60])
    xlabel('time  (min)');
    ylabel('z (px)');

    set(gcf,'position',[1925         123        1558         863],'Color','w')

    %saveas(gcf,strcat(save_path,filesep,char(mice(m)),'_',char(sessions_name(s)),'_pawSessT.png'),'png')
print(gcf, fullfile(save_path, strcat('paw_in_sess_',char(mice(m)),'.pdf')), '-dpdf', '-painters');

end

%% Lick during reach

figure
bw=.2;
axis_xx = [0 2];
sess = 5;
figg=tiledlayout(2,num_animals/2);

for m = 1:num_animals
    nexttile

    plot(tm,squeeze(licks_mm(:,2,:,sess,m)),'color',cat(2,clrs_m(m,:),.5),'linewidth',2); %hold on
    %ylim([0 20])
    set(gca,axeOpt{:})
    xlabel('time (s)'); ylabel('position of tongue (mm)')
    title(mice(m),'color',clrs_m(m,:),'Interpreter','none');
end
title(figg,{'inter-reach interval';char(sessions_name(sess))})
set(gcf,'position',[1950         159        1830         729],'color','w')
%saveas(gcf,strcat(save_path,filesep,'iri',char(sessions_name(sess)),'.png'),'png')


%% Check if there are at least 10 points with tongue detected

win_detection = 5;
is_tongue_detected = squeeze(sum(movsum(~isnan(squeeze(licks_mm(:,1,:,sess,:))),win_detection)==win_detection,1)>1);
trials_lick_detected = sum(is_tongue_detected);
vec_lick = cat(2,trials_lick_detected(wts_range)',trials_lick_detected(a2a_range)');
[h_lick,sig_lick]=ttest2(vec_lick(:,1),vec_lick(:,2));

wd_bx =.5;
   figure
boxplot(vec_lick,'Symbol', 'k.','Color','k','Widths',wd_bx);
hold on
for pos=1:size(vec_lick,2)
    if ismember(pos,1:3:nxh1*3), c=1; else, c=2; end
    f=scatter(ones(1,max_n_group)+(pos-1)+(rand(01,max_n_group)-0.5)/10,vec_lick(:,pos),sz,clrs_m(c,:),'filled','LineWidth',1.5);
end
hold off
%axis square;% ylim([-20 250])
XTickLabel={'CTR';'A2a'};

set(gca,'XTickLabel',XTickLabel); ylabel('nr of trials with licks');
set(gca,axeOpt{:})
title({'Mean end-point in x',sprintf('%s%1.3f','sig = ',sig_lick)})
%set(gcf,'position',[ 2653         452         182         420]);
%print(gcf, fullfile(save_path, 'iri_med_vec.pdf'), '-dpdf', '-painters');

%% TOGETHER

figure
sess = 5;
figg=tiledlayout(2,3);
selectr_ctr = [4,8];
for m = selectr_ctr
    nexttile
    plot(tm,squeeze(licks_mm(:,1,:,sess,m)),'color',cat(2,clrs_m(m,:),.5),'linewidth',2); %hold on
    xlim([tm(1) tm(end)])
    ylim([0 20])
    set(gca,axeOpt{:})
    xlabel('time (s)'); ylabel('position of tongue in x (mm)')
    %title(mice(m),'color',clrs_m(m,:),'Interpreter','none');
end

i=4;
selectr_ctr = [3,9];
for m = selectr_ctr
    nexttile(i)
    plot(tm,squeeze(licks_mm(:,1,:,sess,m)),'color',cat(2,clrs_m(m,:),.5),'linewidth',2); %hold on
    xlim([tm(1) tm(end)])
    ylim([0 20])
    set(gca,axeOpt{:})
    xlabel('time (s)'); ylabel('position of tongue in x(mm)')
    %title(mice(m),'color',clrs_m(m,:),'Interpreter','none');
    i=i+1;
end

nexttile(3,[2 1])
boxplot(vec_lick,'Symbol', 'k.','Color','k','Widths',wd_bx);
hold on
for pos=1:size(vec_lick,2)
    if ismember(pos,1:3:nxh1*3), c=1; else, c=2; end
    f=scatter(ones(1,max_n_group)+(pos-1)+(rand(01,max_n_group)-0.5)/10,vec_lick(:,pos),sz,clrs_m(c,:),'filled','LineWidth',1.5);
end
hold off
%axis square;% ylim([-20 250])
XTickLabel={'CTR';'A2a'};

set(gca,'XTickLabel',XTickLabel); ylabel('nr of trials with licks detected');
set(gca,axeOpt{:})
title({sprintf('%s%1.3f','sig = ',sig_lick)})


set(gcf,'position',[1950         159        1322         729],'color','w')
%saveas(gcf,strcat(save_path,filesep,'iri',char(sessions_name(sess)),'.png'),'png')


print(gcf, fullfile(save_path, 'licks_fig.pdf'), '-dpdf', '-painters');

