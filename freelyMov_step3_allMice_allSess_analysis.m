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


axeOpt = {'linewidth',1.5,'box','off','GridAlpha',0.05,'ticklength',[1,1]*.01};


%%
lw = 1.5;
sessions_TOT = {'S1';'S2';'S3';'S4';'S5'};

figure()
plot(nr_trials_TOT(:,wts_range),'LineWidth',lw,'Color',wts_clr); hold on
plot(nr_trials_TOT(:,a2a_range),'LineWidth',lw,'Color',a2a_clr);hold off
xlim([.8 n_sess_TOT+.2])
set(gca,axeOpt{:})
title('Nr of trials per session')
xticks(1:n_sess_TOT)
xticklabels(sessions_TOT)
xlabel('sessions'); ylabel('number of trials')
legend(cat(1,mice(wts_range),mice(a2a_range)),'box','off','location','southeastoutside','Interpreter', 'none')
saveas(gcf,strcat(save_path,filesep,'number_trials_session_learning.png'),'png')

%% LEARNING
disp_nsess = 3;
figure()
subplot(221)
plot(nr_trials_TOT(:,wts_range),'LineWidth',lw,'Color',wts_clr); hold on
plot(nr_trials_TOT(:,a2a_range),'LineWidth',lw,'Color',a2a_clr);hold off
xlim([.8 n_sess_TOT+.2])
set(gca,axeOpt{:})
title('Nr of trials per session')
xticks(1:n_sess_TOT)
xticklabels(sessions_TOT)
xlabel('sessions'); ylabel('number of trials')
%saveas(gcf,strcat(save_path,filesep,'number_trials_session_learning.png'),'png')

subplot(222)
plot(nr_reaches_t(:,wts_range),'-','LineWidth',lw,'Color',wts_clr); hold on
plot(nr_reaches_t(:,a2a_range),'-','LineWidth',lw,'Color',a2a_clr);hold off

set(gca,axeOpt{:})
title('Nr of reaches per session')
xlim([.9 n_sess+.1])
%ylim([280 1100])
xticks(1:n_sess)
xticklabels(sessions_name)
xlabel('sessions'); ylabel('number of reaches')
%legend(cat(1,mice(wts_range),mice(a2a_range)),'box','off','location','southeastoutside','Interpreter', 'none')
%saveas(gcf,strcat(save_path,filesep,'number_reaches.png'),'png')

subplot(223)
perc_NOTpurpose = squeeze(puorpose_yesNo_nr(2,:,:)./sum(puorpose_yesNo_nr,1));
perc_NOTpurpose(perc_NOTpurpose==1)=nan;

plot(perc_NOTpurpose(:,wts_range),'-o','LineWidth',lw,'Color',wts_clr); hold on
plot(perc_NOTpurpose(:,a2a_range),'-o','LineWidth',lw,'Color',a2a_clr);hold off
xlim([.9 disp_nsess+.1])
set(gca,axeOpt{:})
title('Ratio of non-porpuseful reaches')
xlim([disp_nsess-.1 n_sess+.1])
xticks(disp_nsess:n_sess)
xticklabels(sessions_name(disp_nsess:n_sess))
xlabel('sessions'); ylabel('percentage of non-porpuseful reaches')
legend('','box','off','location','southeastoutside','Interpreter', 'none')

%saveas(gcf,strcat(save_path,filesep,'percent_non_puporse.png'),'png')


subplot(224)
perc_suc = squeeze(suc_yesNo_nr(1,:,:)./puorpose_yesNo_nr(1,:,:));

plot(perc_suc(:,wts_range),'-o','LineWidth',lw,'Color',wts_clr); hold on
plot(perc_suc(:,a2a_range),'-o','LineWidth',lw,'Color',a2a_clr);hold off
set(gca,axeOpt{:})
title('Succcess over purposeful reaches')
xlim([disp_nsess-.1 n_sess+.1])
xticks(disp_nsess:n_sess)
xticklabels(sessions_name(disp_nsess:n_sess))
xlabel('sessions'); ylabel('sucess rate')
legend(cat(1,mice(wts_range),mice(a2a_range)),'box','off','location','southeastoutside','Interpreter', 'none')

%saveas(gcf,strcat(save_path,filesep,'percent_suc_over_purpose.png'),'png')

set(gcf,'Position',[2325         110        1108         832], 'Color','w')
saveas(gcf,strcat(save_path,filesep,'learning.png'),'png')

%%
sess_range = 4:5;
  nxh=1;
  max_n_group = max(length(wts_range),length(a2a_range));

    % Nr suc reaches
    n_suc_reaches_sess = mean(suc_yesNo_nr(1,sess_range,:),2);
    vec_nsuc=nan(max_n_group,(nxh*3)-1);
    vec_nsuc(1:length(wts_range),1:3:nxh*3)=n_suc_reaches_sess(:,wts_range)';
    vec_nsuc(1:length(a2a_range),2:3:nxh*3+1)=n_suc_reaches_sess(:,a2a_range)';
    length_nsuc=size(vec_nsuc,1);
    [h_suc,sig_suc]=ttest2(vec_nsuc(:,1),vec_nsuc(:,2));

figure

nxh=1;
boxplot(vec_nsuc,'Symbol', 'k.','Color','k','Widths',0.8);
hold on
for pos=1:size(vec_nsuc,2)
    if ismember(pos,1:3:nxh*3), c=1; else, c=2; end
    if pos ==1
        f=scatter((1+(rand(01,length_nsuc)-0.5)/8),vec_nsuc(:,pos),80,clrs_m(c,:),'filled','LineWidth',1.5);
    else
        f=scatter(ones(1,length_nsuc)+(pos-1).*(1+(rand(01,length_nsuc)-0.5)/8),vec_nsuc(:,pos),80,clrs_m(c,:),'filled','LineWidth',1.5);
    end
end
hold off
axis square;% ylim([-20 250])
XTickLabel={'CTR';'a2a'};
set(gca,'XTickLabel',XTickLabel,'linewidth',1.5,'box','off','ticklength',[1,1]*.01,'tickdir','out'); 
ylabel('number of successfull reaches');
set(gcf,'color','w')
%title({'Difference','SP: resting onto lifted paw variance', 'and SP: lifted paw variance'})
title({'Nr  reaches',sprintf('%s%1.3f','sig = ',sig_suc)})
saveas(gcf,strcat(save_path,filesep,'nr_suc_reaches.png'),'png')



%% categories

clrs_cats = cat(1,cat1_clr,cat2_clr,cat3_clr);
lw_m = 2;

figure
tiledlayout(2,num_animals/2);

for m=1:num_animals
    nexttile
    plot((cat_reach_nr_t(:,:,m)./sum(cat_reach_nr_t(:,:,m),1))','-o','LineWidth',lw_m)
    xlim([.9 4.1]); ylim([0 1])
    xticks(1:4)
    xticklabels(sessions_name)
    set(gca,axeOpt{:},'colorOrder',clrs_cats)
    title(mice(m),'color',clrs_m(m,:),'Interpreter','none');
end
set(gcf,'position',[1950         159        1830         729],'color','w')
saveas(gcf,strcat(save_path,filesep,'categories_over_session.png'),'png')



%% paw usage R
lw =2;
figure
perc_non_prefPaw = squeeze(sum(r_paw_used_t==0,1))./nr_reaches_t;
plot(perc_non_prefPaw(:,wts_range),'-o','LineWidth',lw,'Color',wts_clr); hold on
plot(perc_non_prefPaw(:,a2a_range),'-o','LineWidth',lw,'Color',a2a_clr);hold off

set(gca,axeOpt{:})
title('Percentage of L paw usage')
xlim([.9 5.1])
xticks(1:4)
xticklabels(sessions_name)
xlabel('sessions'); ylabel('L paw used (%)')
legend(cat(1,mice(wts_range),mice(a2a_range)),'box','off','location','southeastoutside','Interpreter', 'none')

saveas(gcf,strcat(save_path,filesep,'non_dominat_paw.png'),'png')


%% IRI
bw=.2;
axis_xx = [-.5 6];
for sess = 1:num_sessions
    figure
    figg=tiledlayout(2,num_animals/2);
   
    for m = 1:num_animals
        nexttile
        histogram(iri_sec_t(:,sess,m),'binwidth',bw,'Normalization','pdf','facecolor',clrs_m(m,:),'FaceAlpha',0.5); hold on
        xlim(axis_xx)
        set(gca,axeOpt{:})
        xlabel('iri (s)'); ylabel('counts')
        title(mice(m),'color',clrs_m(m,:),'Interpreter','none');
    end
    title(figg,{'inter-reach interval';char(sessions_name(sess))})
    set(gcf,'position',[1950         159        1830         729],'color','w')
    saveas(gcf,strcat(save_path,filesep,'iri',char(sessions_name(sess)),'.png'),'png')

end

%% pooled
figure
figg=tiledlayout(1,1);
sesses = 3:5;
bw=.1;

iri_pool_a2a = iri_sec_t(:,sesses,a2a_range);
iri_pool_wts = iri_sec_t(:,sesses,wts_range);

h1=histogram(iri_pool_wts(~isnan(iri_pool_wts)),'binwidth',bw,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
h2=histogram(iri_pool_a2a(~isnan(iri_pool_a2a)),'binwidth',bw,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold on
xlim(axis_xx)
set(gca,axeOpt{:})

xlabel('iri (s)'); ylabel('pdf')
set(gcf,'color','w','position',[2493         275         603         499])
sess_interval_name = sprintf('%s%i%s%i','S',sesses(1),'-',sesses(end));
title(figg,sprintf('%s%i%s%i','Inter-reach interval S',sesses(1),'-',sesses(end)))

saveas(gcf,strcat(save_path,filesep,'iri_pooled_,',sess_interval_name,'.png'),'png')


%% IRI quantiles

max_quantile = 265;
n_quantiles = 4;

quantiles_idx = nan(max_quantile,n_quantiles,num_sessions,num_animals);
iri_quantiles = nan(max_quantile,n_quantiles,num_sessions,num_animals);

for m=1:num_animals
    for s=1:num_sessions

        iri_tmp = iri_sec_t(:,s,m);
        last_idx = find(isnan(iri_tmp),1);
        quantiles_length = floor(last_idx/n_quantiles);

        for i = 1:n_quantiles
            quantiles_idx(1:quantiles_length,i,s,m) = ((i-1)*quantiles_length)+1 : quantiles_length*i;
            iri_quantiles(1:quantiles_length,i,s,m) = iri_tmp(quantiles_idx(1:quantiles_length,i,s,m));
        end
    end
end

% Mean of quantiles over time
lw = 1.5;
figure

figg=tiledlayout(2,3);
for s=1:num_sessions
nexttile
plot(squeeze(mean(iri_quantiles(:,:,s,wts_range),1,'omitnan')),'-o','linewidth',lw,'color',wts_clr); hold on
plot(squeeze(mean(iri_quantiles(:,:,s,a2a_range),1,'omitnan')),'-o','linewidth',lw,'color',a2a_clr); hold on
set(gca,axeOpt{:})
xticks(1:n_quantiles)
xlabel('quantiles'); ylabel('mean iri (s)')
xlim([.8 4.2])
title(sessions_name(s))
end

set(gcf,'position',[2178         144        1239         837],'color','w');
saveas(gcf,strcat(save_path,filesep,'iri_over_quantiles.png'),'png')


%% IRI across session

iri_sess = cat(1,nan(1,num_sessions,num_animals),diff(reach_times_sess)./fps);
reach_time_min = (reach_times_sess./fps)./60;

m_fit = nan(num_sessions,num_animals);
b_fit = nan(num_sessions,num_animals);

for s=1:num_sessions
    for m=1:num_animals
        x = reach_time_min(2:end,s,m);
        y = iri_sess(2:end,s,m);
        f=fit(x(~isnan(x)),y(~isnan(y)),'poly1');
        m_fit(s,m)=f.p1;
        b_fit(s,m)=f.p2;
    end
end

for s = 1:num_sessions
    figure()
    figg=tiledlayout(2,num_animals/2);
    for m=1:num_animals
        nexttile
        scatter(reach_time_min(:,s,m),iri_sess(:,s,m),8,...
            'filled','MarkerFaceColor',clrs_m(m,:));
        hold on
        plot([0,30],[b_fit(s,m),m_fit(s,m)*30+b_fit(s,m)],...
            'color',clrs_m(m,:),'LineWidth',2)
        set(gca,axeOpt{:})
        xlabel('time in session (min)'); ylabel('inter-reach interval (s)')
        ylim([0 8])
        title({char(mice(m));...
            sprintf('%s%.2f','slope = ',m_fit(s,m))},'Interpreter','none')
    end
    title(figg,{'IRI across session';char(sessions_name(s))})
    set(gcf,'position',[1950         159        1830         729],'color','w')
    saveas(gcf,strcat(save_path,filesep,'iri_across_session_',char(sessions_name(s)),'.png'),'png')
end


%% IRI cumulative across session

for s = 1:num_sessions
    figure()
    %figg=tiledlayout(2,2);
    for m=1:num_animals
        %nexttile

        %max_x = nr_reaches_t(s,m);
        reach_time_min_m = reach_time_min(:,s,m);
        reach_time_min_m_p = reach_time_min_m(~isnan(reach_time_min_m));
        total_reach_nr = length(reach_time_min_m_p);
        idx_reach_norm = (0:total_reach_nr-1)./total_reach_nr;

        plot(idx_reach_norm,reach_time_min_m_p,...
            'color',clrs_m(m,:),'LineWidth',2)
        hold on
        set(gca,axeOpt{:})
        xlabel('reach idx normalized'); ylabel('time in session (min)')
        legend(char(mice),'box','off','location','northwest','Interpreter','none')        %xlim([0 8])
         title(char(sessions_name(s)))
%             sprintf('%s%.2f','slope = ',m_fit(s,m))})
    end
    %title(figg,{'IRI across session';char(sessions_name(s))})
    %set(gcf,'color','w','position',[2212         205         967         752])
    saveas(gcf,strcat(save_path,filesep,'reach_time_evolution_normalized_',char(sessions_name(s)),'.png'),'png')

end

%% IRI pooled
sess_range = 4:5;
x_a2a = reach_time_min(2:end,sess_range,a2a_range);
y_a2a = iri_sess(2:end,sess_range,a2a_range);
f_a2a = fit(x_a2a(~isnan(x_a2a)),y_a2a(~isnan(y_a2a)),'poly1');
m_a2a = f_a2a.p1;
b_a2a = f_a2a.p1;
x_wts = reach_time_min(2:end,sess_range,wts_range);
y_wts = iri_sess(2:end,sess_range,wts_range);
f_wts = fit(x_wts(~isnan(x_wts)),y_wts(~isnan(y_wts)),'poly1');
m_wts = f_wts.p1;
b_wts = f_wts.p1;

figure()
figg=tiledlayout(1,2);

nexttile
scatter(x_a2a(~isnan(x_a2a)),y_a2a(~isnan(y_a2a)),8,...
    'filled','MarkerFaceColor',a2a_clr,'MarkerFaceAlpha',0.5);
hold on
plot([0,30],[b_a2a,m_a2a*30+b_a2a],...
    'color',a2a_clr,'LineWidth',2)
ylim([0 10])
set(gca,axeOpt{:})
xlabel('time in session (min)'); ylabel('inter-reach interval (s)')
            title(sprintf('%s%.2f','A2as slope = ',m_a2a))

nexttile
scatter(x_a2a(~isnan(x_wts)),y_a2a(~isnan(y_wts)),8,...
    'filled','MarkerFaceColor',wts_clr,'MarkerFaceAlpha',0.5);
hold on
plot([0,30],[b_wts,m_wts*30+b_wts],...
    'color',wts_clr,'LineWidth',2)
ylim([0 10])
set(gca,axeOpt{:})
            title(sprintf('%s%.2f','WTs slope = ',m_wts))

xlabel('time in session (min)'); ylabel('inter-reach interval (s)')
title(figg,'IRI across session pooled S6-7')

set(gcf,'Position',[2335         294        1119         502])
saveas(gcf,strcat(save_path,filesep,'iri_across_session_pooled_S6-7.png'),'png')


%% Time from ITI to first reach - peak
bw=.2;
axis_xx = [-.5 6];

for sess = 1:num_sessions
    figure
    figg=tiledlayout(2,num_animals/2);
    for m=1:num_animals
        nexttile
        histogram(time_to_1st_reach_peak_t(:,sess,m),'binwidth',bw,'Normalization','pdf','facecolor',clrs_m(m,:),'FaceAlpha',0.5); hold on
        xlim(axis_xx)
        set(gca,axeOpt{:})
        xlabel('time delay (s)'); ylabel('counts')
    end
    title(figg,{'Time to 1st reach peak after ITI end';char(sessions_name(sess))})
    set(gcf,'position',[1950         159        1830         729],'color','w')

    saveas(gcf,strcat(save_path,filesep,'time_for_1st_reach_peak',char(sessions_name(sess)),'.png'),'png')
end


%% Time from ITI to first reach - start
for sess = 1:num_sessions
    figure
    figg=tiledlayout(2,num_animals/2);
    for m=1:num_animals
        nexttile
        histogram(time_to_1st_reach_start_t(:,sess,m),'binwidth',bw,'Normalization','pdf','facecolor',clrs_m(m,:),'FaceAlpha',0.5); hold on
        xlim(axis_xx)
        set(gca,axeOpt{:})
        xlabel('time delay (s)'); ylabel('pdf')
    end
    title(figg,{'Time to 1st reach start after ITI end';char(sessions_name(sess))})
    set(gcf,'position',[1950         159        1830         729],'color','w')
    %set(gcf,'color','w','WindowStyle', 'Docked')

    saveas(gcf,strcat(save_path,filesep,'time_for_1st_start_peak',char(sessions_name(sess)),'.png'),'png')
end

%% Pooled time from ITI to first reach 
bw=0.1;
axis_xx = [-.5 4];

sess_range = 3:5;
sess_interval_name = sprintf('%s%i%s%i','S',sesses(1),'-',sesses(end));

figure
figg=tiledlayout(1,2);
nexttile
time_to_1st_start_a2a=time_to_1st_reach_start_t(:,sess_range,a2a_range);
time_to_1st_start_wts=time_to_1st_reach_start_t(:,sess_range,wts_range);
histogram(time_to_1st_start_wts(~isnan(time_to_1st_start_wts)),'binwidth',bw,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
histogram(time_to_1st_start_a2a(~isnan(time_to_1st_start_a2a)),'binwidth',bw,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold off
set(gca,axeOpt{:})
xlim(axis_xx)
xlabel('time delay (s)'); ylabel('pdf')
title('time for 1st reach start');

nexttile
time_to_1st_peak_a2a=time_to_1st_reach_peak_t(:,sess_range,a2a_range);
time_to_1st_peak_wts=time_to_1st_reach_peak_t(:,sess_range,wts_range);

histogram(time_to_1st_peak_wts(~isnan(time_to_1st_peak_wts)),'binwidth',bw,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
histogram(time_to_1st_peak_a2a(~isnan(time_to_1st_peak_a2a)),'binwidth',bw,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold off
set(gca,axeOpt{:})
xlim(axis_xx)
xlabel('time delay (s)'); ylabel('pdf')
title('time for 1st reach peak');


legend('wts','a2as','box','off')

set(gcf,'color','w','position',[2420         351        1100         492])
title(figg,strcat('Time to 1st reach after ITI end polled: ',sess_interval_name))
saveas(gcf,strcat(save_path,filesep,'time_for_1st_reach_pooled_S',sess_interval_name,'.png'),'png')


% ----------------------------------------------------------
%% duration of reaches
non_cat3 =  cat_reach_id_t~=3;
dur_reaches_cat12 = dur_reaches_t.*non_cat3;
dur_reaches_cat12(dur_reaches_cat12==0) = nan;
axis_xx = [0 0.35];

for sess = 1:num_sessions
    figure
    figg=tiledlayout(2,num_animals/2);
    for m=1:num_animals
        nexttile
        histogram(dur_reaches_cat12(:,sess,m),'binwidth',0.01,'Normalization','pdf','facecolor',clrs_m(m,:),'FaceAlpha',0.5); hold on
        xlim(axis_xx)
        set(gca,axeOpt{:})
        xlabel('duration forward (s)'); ylabel('counts')
    end
    title(figg,{'Duration of forward reaches';char(sessions_name(sess))})
    set(gcf,'position',[1950         159        1830         729],'color','w')
    pause(1)
    saveas(gcf,strcat(save_path,filesep,'duration_',char(sessions_name(sess)),'.png'),'png')
end

%% duration of reaches over session (by index, not time)
for s = 1:num_sessions
    figure
    figg=tiledlayout(2,num_animals/2);
    for m=1:num_animals
        nexttile
        x_ind = 1:length(dur_reaches_cat12(:,s,m));
        y_dur = dur_reaches_cat12(:,s,m);
        max_x = nr_reaches_t(s,m);

        f=fit(x_ind(~isnan(y_dur))',y_dur(~isnan(y_dur)),'poly1');

        scatter(x_ind,dur_reaches_cat12(:,s,m),8,...
            'filled','MarkerFaceColor',clrs_m(m,:));
        hold on

        plot([0,max_x],[f.p2,f.p1*max_x+f.p2],...
            'color',clrs_m(m,:),'LineWidth',2)
        xlim([0 max_x])
        title({char(mice(m));...
            sprintf('%s%.2f','slope = ',f.p1)},'Interpreter','none')

        %xlim(axis_xx)
        set(gca,axeOpt{:})
        xlabel('reach idx'); ylabel('duration (s)')
    end
    title(figg,{'Duration of forward reaches';char(sessions_name(s))})
    set(gcf,'position',[1950         159        1830         729],'color','w')
    pause(1)
    saveas(gcf,strcat(save_path,filesep,'duration_across_reach',char(sessions_name(sess)),'.png'),'png')
end


%%
range_sess = 4:5;
sess_interval_name = sprintf('%s%i%s%i','S',range_sess(1),'-',range_sess(end));

figure
dur_wts = dur_reaches_cat12(:,range_sess,wts_range);
dur_a2a = dur_reaches_cat12(:,range_sess,a2a_range);
histogram(dur_wts(~isnan(dur_wts)),'binwidth',0.01,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
histogram(dur_a2a(~isnan(dur_a2a)),'binwidth',0.01,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold off
set(gca,axeOpt{:})
xlim(axis_xx)
xlabel('duration forward (s)'); ylabel('counts')
title(sprintf('%s%s','Duration of pooled forward reaches ',sess_interval_name));
legend('wts','a2as','box','off')
saveas(gcf,strcat(save_path,filesep,'durations_pooled_',...
    char(sessions_name(range_sess(1))),'-',char(sessions_name(range_sess(end)))),'png')


%% duration of reaches divided by category
cat1 =  cat_reach_id_t==1;
cat2 =  cat_reach_id_t==2;
dur_reaches_cat1 = dur_reaches_t.*cat1;
dur_reaches_cat2 = dur_reaches_t.*cat2;
dur_reaches_cat1(dur_reaches_cat1==0) = nan;
dur_reaches_cat2(dur_reaches_cat2==0) = nan;
axis_xx = [0 0.35];
bw = 0.01;

for sess = 1:num_sessions
    figure
    figg=tiledlayout(2,num_animals/2);
    for m=1:num_animals
        nexttile
        dur_reaches_cat1_m = dur_reaches_cat1(:,sess,m);
        dur_reaches_cat2_m = dur_reaches_cat2(:,sess,m);
        histogram(dur_reaches_cat2_m(~isnan(dur_reaches_cat2_m)),'binwidth',bw,'Normalization','pdf','facecolor',cat2_clr,'FaceAlpha',0.5); hold on
        histogram(dur_reaches_cat1_m(~isnan(dur_reaches_cat1_m)),'binwidth',bw,'Normalization','pdf','facecolor',cat1_clr,'FaceAlpha',0.5); hold on
        xlim(axis_xx)
        set(gca,axeOpt{:})
        xlabel('duration forward (s)'); ylabel('pdf')
        title(char(mice(m)),'Color',clrs_m(m,:))
    end
    title(figg,{'Duration of forward reaches';char(sessions_name(sess))})
    set(gcf,'position',[1950         159        1830         729],'color','w')
    pause(1)
    saveas(gcf,strcat(save_path,filesep,'duration_byCat_',char(sessions_name(sess)),'.png'),'png')
end

%%
range_sess = 4:5;
sess_interval_name = sprintf('%s%i%s%i','S',range_sess(1),'-',range_sess(end));

dur_wts_c1 = dur_reaches_cat1(:,range_sess,wts_range);
dur_a2a_c1 = dur_reaches_cat1(:,range_sess,a2a_range);
dur_wts_c2 = dur_reaches_cat2(:,range_sess,wts_range);
dur_a2a_c2 = dur_reaches_cat2(:,range_sess,a2a_range);

figure
figg=tiledlayout(2,1);
nexttile
histogram(dur_a2a_c1(~isnan(dur_a2a_c1)),'binwidth',0.01,'Normalization','pdf','facecolor',cat1_clr,'FaceAlpha',0.5); hold on
histogram(dur_a2a_c2(~isnan(dur_a2a_c2)),'binwidth',0.01,'Normalization','pdf','facecolor',cat2_clr,'FaceAlpha',0.5); hold off
set(gca,axeOpt{:})
xlim(axis_xx)
xlabel('duration forward (s)'); ylabel('counts')
title('A2as','color',a2a_clr);
legend('balistic','holding spout','box','off')

nexttile
histogram(dur_wts_c1(~isnan(dur_wts_c1)),'binwidth',0.01,'Normalization','pdf','facecolor',cat1_clr,'FaceAlpha',0.5); hold on
histogram(dur_wts_c2(~isnan(dur_wts_c2)),'binwidth',0.01,'Normalization','pdf','facecolor',cat2_clr,'FaceAlpha',0.5); hold off
set(gca,axeOpt{:})
xlim(axis_xx)
xlabel('duration forward (s)'); ylabel('counts')
title('WTs','color',wts_clr);

    title(figg,'Duration of forward reaches by category')
    set(gcf,'color','w','position',[2554         130         464         840])

saveas(gcf,strcat(save_path,filesep,'durations_buCat_pooled_',...
    char(sessions_name(range_sess(1))),'-',char(sessions_name(range_sess(end)))),'png')
%%

figure
figg=tiledlayout(1,2);
nexttile
histogram(dur_wts_c1(~isnan(dur_wts_c1)),'binwidth',0.01,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
histogram(dur_a2a_c1(~isnan(dur_a2a_c1)),'binwidth',0.01,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold off
set(gca,axeOpt{:})
xlim(axis_xx)
xlabel('duration forward (s)'); ylabel('counts')
title('balistic','color',cat1_clr);
legend('wts','a2as','box','off')

nexttile
histogram(dur_wts_c2(~isnan(dur_wts_c2)),'binwidth',0.01,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
histogram(dur_a2a_c2(~isnan(dur_a2a_c2)),'binwidth',0.01,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold off
set(gca,axeOpt{:})
xlim(axis_xx)
xlabel('duration forward (s)'); ylabel('counts')
title('holding spout','color',cat2_clr);

title(figg,'Duration of forward reaches by category')
set(gcf,'color','w','position',[2273         281        1291         537])

saveas(gcf,strcat(save_path,filesep,'durations_buCat_sep_pooled_',...
    char(sessions_name(range_sess(1))),'-',char(sessions_name(range_sess(end)))),'png')




%% Reaches all SHAPE OVER TIME
transpa = 0.05;
lw_th = 1;
lw_m = 2;
indiv_a_clr = [.8 .8 .8];
x_lim = [0 350];
%z_lim = [100 350];
tm_lim = [tm(1) tm(end)];


for sess = 4:5%1:num_sessions
    figure
    fig2=tiledlayout(4,2);

    for m = 1:4%num_animals
        % m=1;
        nexttile
        plot(tm,squeeze(reaches_all_t(:,1,:,sess,m)),'linewidth',lw_th,'color',cat(2,indiv_a_clr,transpa)); hold on
        plot(tm,median(reaches_all_t(:,1,cat_reach_id_t(:,sess,m)==3,sess,m),3,'omitnan'),'linewidth',lw_m,'color',cat3_clr); hold on
        plot(tm,median(reaches_all_t(:,1,cat_reach_id_t(:,sess,m)==2,sess,m),3,'omitnan'),'linewidth',lw_m,'color',cat2_clr); hold on
        plot(tm,median(reaches_all_t(:,1,cat_reach_id_t(:,sess,m)==1,sess,m),3,'omitnan'),'linewidth',lw_m,'color',cat1_clr); hold on
        xlim(tm_lim); ylim(x_lim);
        xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
        set(gca,axeOpt{:})
        xlabel('time (s)'); ylabel('x (px)')
        title(mice(m),'Color',clrs_m(m,:))

        if rem(m,2)~=0, nexttile(m*2+1), else, nexttile, end
        plot(tm,squeeze(reaches_all_t(:,2,:,sess,m)),'linewidth',lw_th,'color',cat(2,indiv_a_clr,transpa)); hold on
        plot(tm,median(reaches_all_t(:,2,cat_reach_id_t(:,sess,m)==3,sess,m),3,'omitnan'),'linewidth',lw_m,'color',cat3_clr); hold on
        plot(tm,median(reaches_all_t(:,2,cat_reach_id_t(:,sess,m)==2,sess,m),3,'omitnan'),'linewidth',lw_m,'color',cat2_clr); hold on
        plot(tm,median(reaches_all_t(:,2,cat_reach_id_t(:,sess,m)==1,sess,m),3,'omitnan'),'linewidth',lw_m,'color',cat1_clr); hold on
        xlim(tm_lim); ylim(x_lim);
        xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
        set(gca,axeOpt{:})
        xlabel('time (s)'); ylabel('z (px)')

    end
    title(fig2,{'paw position over time';char(sessions_name(sess))})
    set(gcf,'position',[ 2463         147        1006         848],'Color','w');

    saveas(gcf,strcat(save_path,filesep,'xz_overtime_',char(sessions_name(sess)),'.png'),'png')

end


%%  DISTRIBUTION OF X AND Z POSITION during forward
max_reach = 41;
range_histogram = 1:max_reach;
bw = 10;
axis_xx = [50 450];

for sess = 1:num_sessions
    figure
    figg=tiledlayout(2,num_animals/2);

    for m = 1:num_animals
        nexttile
        histogram(reaches_all_t(range_histogram,2,:,sess,m),'binwidth',bw,'Normalization','pdf','facecolor',clrs_m(m,:),'FaceAlpha',0.5); hold on
        xlim(axis_xx)
        set(gca,axeOpt{:})
        xlabel('x (px)'); ylabel('counts')
        title(char(mice(m)),'Color',clrs_m(m,:),'Interpreter','none')

    end
    title(figg,{'z position';char(sessions_name(sess))})
    set(gcf,'position',[1950         159        1830         729],'color','w')

    saveas(gcf,strcat(save_path,filesep,'distribution_of_position_z',char(sessions_name(sess)),'.png'),'png')

end

%%
figure
figg=tiledlayout(1,2);
bw = 10;
axis_xx = [0 300];
axis_zz = [50 450];
sesses = 4:5;
sess_interval_name = sprintf('%s%i%s%i','S',sesses(1),'-',sesses(end));

nexttile
pos_x_wts = reaches_all_t(range_histogram,1,:,sesses,wts_range);
pos_x_a2a = reaches_all_t(range_histogram,1,:,sesses,a2a_range);
h1=histogram(pos_x_wts(~isnan(pos_x_wts)),'binwidth',bw,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
h2=histogram(pos_x_a2a(~isnan(pos_x_a2a)),'binwidth',bw,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold on
xlim(axis_xx)
set(gca,axeOpt{:})
xlabel('x (px)'); ylabel('counts')
legend('wts','a2as','box','off')

nexttile
pos_z_wts = reaches_all_t(range_histogram,2,:,sesses,wts_range);
pos_z_a2a = reaches_all_t(range_histogram,2,:,sesses,a2a_range);
h1=histogram(pos_z_wts(~isnan(pos_z_wts)),'binwidth',bw,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
h2=histogram(pos_z_a2a(~isnan(pos_z_a2a)),'binwidth',bw,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold on
xlim(axis_zz)
set(gca,axeOpt{:})
xlabel('z (px)'); ylabel('pdf')

set(gcf,'color','w','position',[2646         351         874         368])

title(figg,sprintf('%s%s','Distribution of paw position from time 0 to reach ',sess_interval_name))
saveas(gcf,strcat(save_path,filesep,'distribution_of_position_XZ_pooled_',sess_interval_name,'.png'),'png')

%% END-POINT IN X AND Z
figure
figg=tiledlayout(1,2);
bw = 5;
%axis_xx = [190 350];
range_histogram = 1:61;
sesses = 3:5;
sess_interval_name = sprintf('%s%i%s%i','S',sesses(1),'-',sesses(end));

nexttile
end_point_x_wt = max(reaches_all_t(range_histogram,1,:,sesses,wts_range));
end_point_x_a2 = max(reaches_all_t(range_histogram,1,:,sesses,a2a_range));
h1=histogram(end_point_x_wt(~isnan(end_point_x_wt)),'binwidth',bw,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
h2=histogram(end_point_x_a2(~isnan(end_point_x_a2)),'binwidth',bw,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold on
%xlim(axis_xx)
set(gca,axeOpt{:})
xlabel('end point in x (px)'); ylabel('pdf')
legend('wts','a2as','box','off')

nexttile
bw = 10;
end_point_z_wt = max(reaches_all_t(range_histogram,2,:,sesses,wts_range));
end_point_z_a2 = max(reaches_all_t(range_histogram,2,:,sesses,a2a_range));
h1=histogram(end_point_z_wt(~isnan(end_point_z_wt)),'binwidth',bw,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
h2=histogram(end_point_z_a2(~isnan(end_point_z_a2)),'binwidth',bw,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold on
xlim([150 450])
set(gca,axeOpt{:})
xlabel('end point in z (px)'); ylabel('pdf')
legend('wts','a2as','box','off')

set(gcf,'color','w','position',[2646         351         874         368])
title(figg,sprintf('%s%s%s%s','Distribution of end points pooled ' ,...
    char(sessions_name(sesses(1))),'-',char(sessions_name(sesses(end)))))

saveas(gcf,strcat(save_path,filesep,'endpoints_pooled_',...
    char(sessions_name(sesses(1))),'-',char(sessions_name(sesses(end)))),'png')


%% Reaches all SPEED OVER TIME
transpa = 0.05;
lw_th = 1;
lw_m = 2;
indiv_a_clr = [.8 .8 .8];
indiv_w_clr = [.8 .8 .8];
x_lim = [-20 30];
z_lim = [-20 30];
tm_lim = [tm(1) tm(end)];


for sess = 4:5%1:num_sessions

    figure
    fig2=tiledlayout(4,2);

    for m = 1:4%num_animals
        nexttile
        plot(tm,squeeze(speed_reaches_all_t(:,1,:,sess,m)),'linewidth',lw_th,'color',cat(2,indiv_a_clr,transpa)); hold on
        plot(tm,median(speed_reaches_all_t(:,1,cat_reach_id_t(:,sess,m)==3,sess,m),3,'omitnan'),'linewidth',lw_m,'color',cat3_clr); hold on
        plot(tm,median(speed_reaches_all_t(:,1,cat_reach_id_t(:,sess,m)==2,sess,m),3,'omitnan'),'linewidth',lw_m,'color',cat2_clr); hold on
        plot(tm,median(speed_reaches_all_t(:,1,cat_reach_id_t(:,sess,m)==1,sess,m),3,'omitnan'),'linewidth',lw_m,'color',cat1_clr); hold on
        xlim(tm_lim); ylim(x_lim);
        xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
        set(gca,axeOpt{:})
        xlabel('time (s)'); ylabel('x (px)')
        title(mice(m),'Color',clrs_m(m,:))

        if rem(m,2)~=0, nexttile(m*2+1), else, nexttile, end
        plot(tm,squeeze(speed_reaches_all_t(:,2,:,sess,m)),'linewidth',lw_th,'color',cat(2,indiv_a_clr,transpa)); hold on
        plot(tm,median(speed_reaches_all_t(:,2,cat_reach_id_t(:,sess,m)==3,sess,m),3,'omitnan'),'linewidth',lw_m,'color',cat3_clr); hold on
        plot(tm,median(speed_reaches_all_t(:,2,cat_reach_id_t(:,sess,m)==2,sess,m),3,'omitnan'),'linewidth',lw_m,'color',cat2_clr); hold on
        plot(tm,median(speed_reaches_all_t(:,2,cat_reach_id_t(:,sess,m)==1,sess,m),3,'omitnan'),'linewidth',lw_m,'color',cat1_clr); hold on
        xlim(tm_lim); ylim(x_lim);
        xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
        set(gca,axeOpt{:})
        xlabel('time (s)'); ylabel('z (px)')
    end

    title(fig2,{'paw speed over time';char(sessions_name(sess))})
    set(gcf,'position',[ 2463         147        1006         848],'Color','w');

    saveas(gcf,strcat(save_path,filesep,'speed_xz_overtime_',char(sessions_name(sess)),'.png'),'png')
end


%% Reaches all DISTANCE & 2D SPEED
transpa = 0.05;
lw_th = 1;
lw_m = 2;
indiv_a_clr = [.8 .8 .8];
%indiv_w_clr = [.8 .8 .8];
x_lim = [0 250];
z_lim = [-20 30];
tm_lim = [tm(1) tm(end)];


for sess = 1:num_sessions

    figure
    fig2=tiledlayout(4,2);

    for m = 1:num_animals
        nexttile
        plot(tm,squeeze(distance_water_all_t(:,:,sess,m)),'linewidth',lw_th,'color',cat(2,indiv_a_clr,transpa)); hold on
        plot(tm,median(distance_water_all_t(:,cat_reach_id_t(:,sess,m)==3,sess,m),2,'omitnan'),'linewidth',lw_m,'color',cat3_clr); hold on
        plot(tm,median(distance_water_all_t(:,cat_reach_id_t(:,sess,m)==2,sess,m),2,'omitnan'),'linewidth',lw_m,'color',cat2_clr); hold on
        plot(tm,median(distance_water_all_t(:,cat_reach_id_t(:,sess,m)==1,sess,m),2,'omitnan'),'linewidth',lw_m,'color',cat1_clr); hold on
        xlim(tm_lim); ylim(x_lim);
        xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
        set(gca,axeOpt{:})
        xlabel('time (s)'); ylabel('distance to water (px)')
        title(mice(m),'Color',clrs_m(m,:))

        if rem(m,2)~=0, nexttile(m*2+1), else, nexttile, end
        plot(tm,squeeze(speed_water_all_t(:,:,sess,m)),'linewidth',lw_th,'color',cat(2,indiv_a_clr,transpa)); hold on
        plot(tm,median(speed_water_all_t(:,cat_reach_id_t(:,sess,m)==3,sess,m),2,'omitnan'),'linewidth',lw_m,'color',cat3_clr); hold on
        plot(tm,median(speed_water_all_t(:,cat_reach_id_t(:,sess,m)==2,sess,m),2,'omitnan'),'linewidth',lw_m,'color',cat2_clr); hold on
        plot(tm,median(speed_water_all_t(:,cat_reach_id_t(:,sess,m)==1,sess,m),2,'omitnan'),'linewidth',lw_m,'color',cat1_clr); hold on
        xlim(tm_lim); ylim(z_lim);
        set(gca,axeOpt{:})
        xlabel('time (s)'); ylabel('speed towards water (px/s)')
        xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
        set(gca,axeOpt{:})

    end

    title(fig2,{'Distance and 2D Speed to water';char(sessions_name(sess))})
    set(gcf,'position',[ 2463         147        1006         848],'Color','w');

    saveas(gcf,strcat(save_path,filesep,'distance_and_2dspeed',char(sessions_name(sess)),'.png'),'png')

end

%%
lw_m = 1;
figure
plot(tm,squeeze(median(speed_water_all_t(:,cat_reach_id_t(:,sess,m)==1,sess,wts_range),2,'omitnan')),'linewidth',lw_m,'color',wts_clr); hold on
plot(tm,squeeze(median(speed_water_all_t(:,cat_reach_id_t(:,sess,m)==1,sess,a2a_range),2,'omitnan')),'linewidth',lw_m,'color',a2a_clr); hold on

plot(tm,mean(squeeze(median(speed_water_all_t(:,cat_reach_id_t(:,sess,m)==1,sess,wts_range),2,'omitnan')),2),'linewidth',2,'color',wts_clr); hold on
plot(tm,mean(squeeze(median(speed_water_all_t(:,cat_reach_id_t(:,sess,m)==1,sess,a2a_range),2,'omitnan')),2),'linewidth',2,'color',a2a_clr); hold on
xlabel('time (s)'); ylabel('speed towards water (px/s)')
xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
set(gca,axeOpt{:})
set(gcf,'position',[2488         362         968         407],'Color','w');

saveas(gcf,strcat(save_path,filesep,'2dspeed_all',char(sessions_name(sess)),'.png'),'png')



%%
max_reach = 41;
range_histogram = 1:max_reach;
bw = 2;
axis_xx = [-35 35];

for sess = 1:num_sessions
    figure
    figg=tiledlayout(2,num_animals/2);

    
    for m = 1:num_animals
        nexttile
        histogram(speed_reaches_all_t(range_histogram,1,:,sess,m),'binwidth',bw,'Normalization','pdf','facecolor',clrs_m(m,:),'FaceAlpha',0.5); hold on
        xlim(axis_xx)
        set(gca,axeOpt{:})
        xlabel('speed in x (px/s)'); ylabel('pdf')
    end

    title(figg,{'partial velocity in x - reach to peak';char(sessions_name(sess))})
    set(gcf,'position',[1950         159        1830         729],'color','w')

    saveas(gcf,strcat(save_path,filesep,'distribution_of_speed_x_untilReach',char(sessions_name(sess)),'.png'),'png')

end


%%
bw = 2;
axis_xx = [-35 35];

for sess = 1:num_sessions
    figure
    figg=tiledlayout(2,num_animals/2);

    for m=1:num_animals
    nexttile
    histogram(speed_reaches_all_t(range_histogram,2,:,sess,m),'binwidth',bw,'Normalization','pdf','facecolor',clrs_m(m,:),'FaceAlpha',0.5); hold on
    xlim(axis_xx)
    set(gca,axeOpt{:})
    xlabel('speed in z (px/s)'); ylabel('pdf')
    end
  

    title(figg,{'partial velocity in z - reach to peak';char(sessions_name(sess))})
    set(gcf,'position',[1950         159        1830         729],'color','w')
    saveas(gcf,strcat(save_path,filesep,'distribution_of_speed_z_untilReach',char(sessions_name(sess)),'.png'),'png')



end

%%
figure
figg=tiledlayout(1,2);

sesses = 4:5;

nexttile
speed_x_a2a = speed_reaches_all_t(range_histogram,1,:,sesses,a2a_range);
speed_x_wts = speed_reaches_all_t(range_histogram,1,:,sesses,wts_range);
histogram(speed_x_wts(~isnan(speed_x_wts)),'binwidth',bw,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
histogram(speed_x_a2a(~isnan(speed_x_a2a)),'binwidth',bw,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold on
xlim(axis_xx)
set(gca,axeOpt{:})
xlabel('partial velocity x (px/s)'); ylabel('pdf')

nexttile
speed_z_a2a = speed_reaches_all_t(range_histogram,2,:,sesses,a2a_range);
speed_z_wts = speed_reaches_all_t(range_histogram,2,:,sesses,wts_range);
histogram(speed_z_wts(~isnan(speed_z_wts)),'binwidth',bw,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
histogram(speed_z_a2a(~isnan(speed_z_a2a)),'binwidth',bw,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold on
xlim(axis_xx)
set(gca,axeOpt{:})
xlabel('partial velocity z (px/s)'); ylabel('pdf')
legend('wts','a2as','box','off')

set(gcf,'color','w','position',[2646         351         874         368])
title(figg,'Distribution of partial velocities from time 0 to reach (S4-5)')
saveas(gcf,strcat(save_path,filesep,'distribution_of_speed_pooled_S4-5.png'),'png')


%% HISTOGRAM WITH JUST REACH FORWARD
speed_r_f_only = nan(41,2,n_sess,num_animals);

for m=1:num_animals
    for s=1:n_sess
        for rr = 1:nr_reaches_max
            if ~isnan(start_fwd_t(rr,s,m))
                speed_r_f_only(start_fwd_t(rr,s,m):41,:,rr,s,m) = speed_reaches_all_t(start_fwd_t(rr,s,m):41,:,rr,s,m);
            else
                speed_r_f_only(:,:,rr,s,m) = nan(41,2);
            end

        end
    end
end
speed_r_f_only(speed_r_f_only==0)=nan;

%%

for sess = 1:num_sessions
    figure
    figg=tiledlayout(2,num_animals/2);

    
for m=1:num_animals
    nexttile
    histogram(speed_r_f_only(:,1,:,sess,m),'binwidth',bw,'Normalization','pdf','facecolor',clrs_m(m,:),'FaceAlpha',0.5); hold on
    xlim(axis_xx)
    set(gca,axeOpt{:})
    xlabel('partial velocity x (px/s)'); ylabel('counts')
end


    title(figg,{'partial velocity on x of forward movement only during reach';char(sessions_name(sess))})
    set(gcf,'position',[1950         159        1830         729],'color','w')
    saveas(gcf,strcat(save_path,filesep,'distribution_of_speed_forward_x_onlyREACH',char(sessions_name(sess)),'.png'),'png')

end


%%

figure
figg=tiledlayout(1,2);
bw = 2;
axis_xx_s = [-20 35];
axis_zz_s = [-30 40];
sesses = 4:5;

nexttile
speed_x_f_a2a = speed_r_f_only(:,1,:,sesses,a2a_range);
speed_x_f_wts = speed_r_f_only(:,1,:,sesses,wts_range);
histogram(speed_x_f_wts(~isnan(speed_x_f_wts)),'binwidth',bw,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
histogram(speed_x_f_a2a(~isnan(speed_x_f_a2a)),'binwidth',bw,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold on
xlim(axis_xx_s)
set(gca,axeOpt{:})
xlabel('velocity in x (px/s)'); ylabel('counts')
legend('wts','a2as','box','off')

nexttile
speed_z_f_a2a = speed_r_f_only(:,2,:,sesses,a2a_range);
speed_z_f_wts = speed_r_f_only(:,2,:,sesses,wts_range);
histogram(speed_z_f_wts(~isnan(speed_z_f_wts)),'binwidth',bw,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
histogram(speed_z_f_a2a(~isnan(speed_z_f_a2a)),'binwidth',bw,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold on
xlim(axis_zz_s)
set(gca,axeOpt{:})
xlabel('velocity in z (px/s)'); ylabel('counts')

set(gcf,'color','w','position',[2646         351         874         368])
title(figg,'Distribution of partial velocity only during the reach (S4-5)')

saveas(gcf,strcat(save_path,filesep,'distribution_of_speed_onlyREACH_pooled_S4-5',char(sessions_name(sess)),'.png'),'png')





%% Trajectories


transpa = 0.05;
lw_th = 1;
lw_m = 2;
indiv_a_clr = [.8 .8 .8];
indiv_w_clr = [.8 .8 .8];
x_lim = [50 320];
z_lim = [80 350];
tm_lim = [tm(1) tm(end)];


for sess = 4:5%1:num_sessions

    figure
    fig2=tiledlayout(2,2);

for m=1:4%num_animals
    nexttile
    plot(squeeze(reaches_all_t(:,1,:,sess,m)),squeeze(reaches_all_t(:,2,:,sess,m)),'linewidth',lw_th,'color',cat(2,indiv_a_clr,transpa)); hold on
    plot(median(reaches_all_t(:,1,cat_reach_id_t(:,sess,m)==3,sess,m),3,'omitnan'),median(reaches_all_t(:,2,cat_reach_id_t(:,sess,m)==3,sess,m),3,'omitnan'),'linewidth',lw_m,'color',cat3_clr); hold on
    plot(median(reaches_all_t(:,1,cat_reach_id_t(:,sess,m)==2,sess,m),3,'omitnan'),median(reaches_all_t(:,2,cat_reach_id_t(:,sess,m)==2,sess,m),3,'omitnan'),'linewidth',lw_m,'color',cat2_clr); hold on
    plot(median(reaches_all_t(:,1,cat_reach_id_t(:,sess,m)==1,sess,m),3,'omitnan'),median(reaches_all_t(:,2,cat_reach_id_t(:,sess,m)==1,sess,m),3,'omitnan'),'linewidth',lw_m,'color',cat1_clr); hold on

    plot(median(reaches_all_t(41,1,cat_reach_id_t(:,sess,m)==3,sess,m),3,'omitnan'),median(reaches_all_t(41,2,cat_reach_id_t(:,sess,m)==3,sess,m),3,'omitnan'),'o','linewidth',5,'color',cat3_clr); hold on
    plot(median(reaches_all_t(41,1,cat_reach_id_t(:,sess,m)==2,sess,m),3,'omitnan'),median(reaches_all_t(41,2,cat_reach_id_t(:,sess,m)==2,sess,m),3,'omitnan'),'o','linewidth',5,'color',cat2_clr); hold on
    plot(median(reaches_all_t(41,1,cat_reach_id_t(:,sess,m)==1,sess,m),3,'omitnan'),median(reaches_all_t(41,2,cat_reach_id_t(:,sess,m)==1,sess,m),3,'omitnan'),'o','linewidth',5,'color',cat1_clr); hold on

    xlim(x_lim); ylim(z_lim)
    xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
    set(gca,axeOpt{:})
    xlabel('x (px)'); ylabel('z (px)')
    title(mice(m),'Color',clrs_m(m,:))
end

    title(fig2,{'paw position over time';char(sessions_name(sess))})
    set(gcf,'position',[ 2100          84        1080         872],'Color','w');


    saveas(gcf,strcat(save_path,filesep,'xz_trajectories_full',char(sessions_name(sess)),'.png'),'png')

end


%% Trajectories

range_tj = 1:41;
transpa = 0.05;
lw_th = 1;
lw_m = 2;
indiv_a_clr = [.8 .8 .8];
indiv_w_clr = [.8 .8 .8];
x_lim = [50 320];
z_lim = [80 350];
tm_lim = [tm(1) tm(end)];


for sess = 4:5%1:num_sessions

    figure
    fig2=tiledlayout(2,2);

    for m=1:4%num_animals
        nexttile
        plot(squeeze(reaches_all_t(range_tj,1,:,sess,m)),squeeze(reaches_all_t(range_tj,2,:,sess,m)),'linewidth',lw_th,'color',cat(2,indiv_a_clr,transpa)); hold on
        plot(median(reaches_all_t(range_tj,1,cat_reach_id_t(:,sess,m)==3,sess,m),3,'omitnan'),median(reaches_all_t(range_tj,2,cat_reach_id_t(:,sess,m)==3,sess,m),3,'omitnan'),'linewidth',lw_m,'color',cat3_clr); hold on
        plot(median(reaches_all_t(range_tj,1,cat_reach_id_t(:,sess,m)==2,sess,m),3,'omitnan'),median(reaches_all_t(range_tj,2,cat_reach_id_t(:,sess,m)==2,sess,m),3,'omitnan'),'linewidth',lw_m,'color',cat2_clr); hold on
        plot(median(reaches_all_t(range_tj,1,cat_reach_id_t(:,sess,m)==1,sess,m),3,'omitnan'),median(reaches_all_t(range_tj,2,cat_reach_id_t(:,sess,m)==1,sess,m),3,'omitnan'),'linewidth',lw_m,'color',cat1_clr); hold on

        plot(median(reaches_all_t(41,1,cat_reach_id_t(:,sess,m)==3,sess,m),3,'omitnan'),median(reaches_all_t(41,2,cat_reach_id_t(:,sess,m)==3,sess,m),3,'omitnan'),'o','linewidth',5,'color',cat3_clr); hold on
        plot(median(reaches_all_t(41,1,cat_reach_id_t(:,sess,m)==2,sess,m),3,'omitnan'),median(reaches_all_t(41,2,cat_reach_id_t(:,sess,m)==2,sess,m),3,'omitnan'),'o','linewidth',5,'color',cat2_clr); hold on
        plot(median(reaches_all_t(41,1,cat_reach_id_t(:,sess,m)==1,sess,m),3,'omitnan'),median(reaches_all_t(41,2,cat_reach_id_t(:,sess,m)==1,sess,m),3,'omitnan'),'o','linewidth',5,'color',cat1_clr); hold on

        xlim(x_lim); ylim(z_lim)
        xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
        set(gca,axeOpt{:})
        xlabel('x (px)'); ylabel('z (px)')
        title(mice(m),'Color',clrs_m(m,:))
    end
    title(fig2,{'paw position over time';char(sessions_name(sess))})
    set(gcf,'position',[ 2100          84        1080         872],'Color','w');

    saveas(gcf,strcat(save_path,filesep,'xz_trajectories_until_reach',char(sessions_name(sess)),'.png'),'png')

end



%% Check distance to mean

cat1_y =  cat_reach_id_t==1;
cat1_all = reaches_all_t.*permute(repmat(cat1_y,[1,1,1,2,101]),[5,4,1,2,3]);
cat1_all(cat1_all==0) = nan;

cat2_y =  cat_reach_id_t==2;
cat2_all = reaches_all_t.*permute(repmat(cat2_y,[1,1,1,2,101]),[5,4,1,2,3]);
cat2_all(cat2_all==0) = nan;

cat1_all_med = squeeze(median(cat1_all,3,'omitnan'));
cat2_all_med = squeeze(median(cat2_all,3,'omitnan'));
%axis_xx = [0 0.35];


% SPACE-ONLY -------------------------

n_timepoints = size(reaches_all_t,1);
% Category 1
cat1_all_distMed_sp = nan(n_timepoints,nr_reaches_max,num_sessions,num_animals);

for m = 1: num_animals
    for s = 1:num_sessions
        for tt=1:nr_reaches_max
            for pt=1:n_timepoints
                % TO  THE MEAN
                tmp_dist_allMean=sqrt((cat1_all(pt,1,tt,s,m)-cat1_all_med(:,1,s,m)).^2+...
                    (cat1_all(pt,2,tt,s,m)-cat1_all_med(:,2,s,m)).^2);

                [min_dist_mean,min_dist_loc_mean]=min(tmp_dist_allMean);
                cat1_all_distMed_sp(pt,tt,s,m)=min_dist_mean;
            end
        end
    end
end

% Category 2
cat2_all_distMed_sp = nan(n_timepoints,nr_reaches_max,num_sessions,num_animals);

for m = 1: num_animals
    for s = 1:num_sessions
        for tt=1:nr_reaches_max
            for pt=1:n_timepoints
                % TO  THE MEAN
                tmp_dist_allMean=sqrt((cat2_all(pt,1,tt,s,m)-cat2_all_med(:,1,s,m)).^2+...
                    (cat2_all(pt,2,tt,s,m)-cat2_all_med(:,2,s,m)).^2);

                [min_dist_mean,min_dist_loc_mean]=min(tmp_dist_allMean);
                cat2_all_distMed_sp(pt,tt,s,m)=min_dist_mean;
            end
        end
    end
end

%% HISTOGRAM

bw = 2;
axis_xx = [-5, 80];
range_hist = 20:41;

for sess = 1:num_sessions
    figure
    figg=tiledlayout(2,num_animals/2);

    for m=1:num_animals
    nexttile
    histogram(cat1_all_distMed_sp(range_hist,:,sess,m),'binwidth',bw,'Normalization','pdf','facecolor',clrs_m(m,:),'FaceAlpha',0.5); hold on
    xlim(axis_xx)
    set(gca,axeOpt{:})
    xlabel('distance to med (px)'); ylabel('counts')
    end

    title(figg,{'integram of distance to median cat 1';char(sessions_name(sess))})
    set(gcf,'position',[1950         159        1830         729],'color','w')
    saveas(gcf,strcat(save_path,filesep,'distance_to_MED_cat1_aroundReach',char(sessions_name(sess)),'.png'),'png')

end



%% HISTOGRAM cat2

bw = 2;
axis_xx = [-5, 80];
range_hist = 20:41;

for sess = 1:num_sessions
    figure
    figg=tiledlayout(2,num_animals/2);

    for m=1:num_animals
    nexttile
    histogram(cat2_all_distMed_sp(range_hist,:,sess,m),'binwidth',bw,'Normalization','pdf','facecolor',clrs_m(m,:),'FaceAlpha',0.5); hold on
    xlim(axis_xx)
    set(gca,axeOpt{:})
    xlabel('distance to med (px)'); ylabel('counts')
    end

    title(figg,{'integram of distance to median cat 2';char(sessions_name(sess))})
    set(gcf,'position',[1950         159        1830         729],'color','w')
    saveas(gcf,strcat(save_path,filesep,'distance_to_MED_cat2_aroundReach',char(sessions_name(sess)),'.png'),'png')

end

%%
figure
figg=tiledlayout(1,2);

sesses = 3:4;

nexttile
range_hist1 = 1:101;
cat1_disMed_full_a2a = cat1_all_distMed_sp(range_hist1,:,sesses,a2a_range);
cat1_disMed_full_wts = cat1_all_distMed_sp(range_hist1,:,sesses,wts_range);
histogram(cat1_disMed_full_wts(~isnan(cat1_disMed_full_wts)),'binwidth',bw,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
histogram(cat1_disMed_full_a2a(~isnan(cat1_disMed_full_a2a)),'binwidth',bw,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold on
xlim(axis_xx)
set(gca,axeOpt{:})
xlabel('distance to median (px)'); ylabel('pdf')
title('full range')

nexttile
range_hist2 = 30:41;
cat1_disMed_r_a2a = cat1_all_distMed_sp(range_hist2,:,sesses,a2a_range);
cat1_disMed_r_wts = cat1_all_distMed_sp(range_hist2,:,sesses,wts_range);
histogram(cat1_disMed_r_wts(~isnan(cat1_disMed_r_wts)),'binwidth',bw,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
histogram(cat1_disMed_r_a2a(~isnan(cat1_disMed_r_a2a)),'binwidth',bw,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold on
xlim(axis_xx)
set(gca,axeOpt{:})
xlabel('distance to median (px)'); ylabel('pdf')
title('range around reach')
legend('wts','a2as','box','off')

set(gcf,'color','w','position',[2646         351         874         368])
title(figg, sprintf('%s%s%s%s',...
    'distance to median of balistic reaches ',...
    char(sessions_name(sesses(1))),'-',char(sessions_name(sesses(end)))))

saveas(gcf,strcat(save_path,filesep,'distribution_of_distance_to_median_pooled_cat1',...
    char(sessions_name(sesses(1))),'-',char(sessions_name(sesses(end)))),'png')





%%
figure
figg=tiledlayout(1,2);

sesses = 3:4;

nexttile
range_hist1 = 1:101;
cat2_disMed_full_a2a = cat2_all_distMed_sp(range_hist1,:,sesses,a2a_range);
cat2_disMed_full_wts = cat2_all_distMed_sp(range_hist1,:,sesses,wts_range);
histogram(cat2_disMed_full_wts(~isnan(cat2_disMed_full_wts)),'binwidth',bw,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
histogram(cat2_disMed_full_a2a(~isnan(cat2_disMed_full_a2a)),'binwidth',bw,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold on
xlim(axis_xx)
set(gca,axeOpt{:})
xlabel('distance to median (px)'); ylabel('pdf')
title('full range')

nexttile
range_hist2 = 30:41;
cat2_disMed_r_a2a = cat2_all_distMed_sp(range_hist2,:,sesses,a2a_range);
cat2_disMed_r_wts = cat2_all_distMed_sp(range_hist2,:,sesses,wts_range);
histogram(cat2_disMed_r_wts(~isnan(cat2_disMed_r_wts)),'binwidth',bw,'Normalization','pdf','facecolor',wts_clr,'FaceAlpha',0.5); hold on
histogram(cat2_disMed_r_a2a(~isnan(cat2_disMed_r_a2a)),'binwidth',bw,'Normalization','pdf','facecolor',a2a_clr,'FaceAlpha',0.5); hold on
xlim(axis_xx)
set(gca,axeOpt{:})
xlabel('distance to median (px)'); ylabel('pdf')
title('range around reach')
legend('wts','a2as','box','off')

set(gcf,'color','w','position',[2646         351         874         368])
title(figg, sprintf('%s%s%s%s',...
    'distance to median of holding spout reaches ',...
    char(sessions_name(sesses(1))),'-',char(sessions_name(sesses(end)))))

saveas(gcf,strcat(save_path,filesep,'distribution_of_distance_to_median_pooled_cat2',...
    char(sessions_name(sesses(1))),'-',char(sessions_name(sesses(end)))),'png')






