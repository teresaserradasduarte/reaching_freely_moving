%% Session Analysis STEP2
% After extracting an categorizing reaches, display the data
% v3 -> v4 : add/replace origin from water to centre of resting bar
clear; close all; clc

%% Load side view
raw_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\raw_data';
mat_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\analyzed_data\mat_files';
out_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\analyzed_data\output_files';

group = '20250106_A2aCasp_G2';
setup = 'freely_mov';

save_flag = true;

% Allocate space for mouse cells
data_path = (strcat(mat_folder,filesep,group,filesep,setup));
folders_mice = dir(data_path);
folders_mice=folders_mice(3:end,:);
num_animals=size(folders_mice,1);

% Allocate space for mouse cells
mice = cell(num_animals,1);
mice_path = cell(num_animals,1);
%num_sessions_tmp = 4;

for m = 7:10%num_animals
%m=1;
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
   %s=5;
        sessions_name(s,1)=cellstr(convertCharsToStrings(folders_session(s,1).name));
        sessions_path(s,1) = strcat(data_path, filesep, mice(m,1), filesep, sessions_name(s,1), filesep);

         fprintf('%s%s%s%s%s','loading mouse ',char(mice(m,1)),', session ',char(sessions_name(s,1)))
            fprintf('\n')

        mouse = char(mice(m,1));
        session = char(sessions_name(s,1));

        %% Load data of single session
        load(strcat(mat_folder,filesep,group,filesep,setup,filesep,mouse,filesep,session,filesep,'session_reaching_data.mat')) % mouse_info
        sess_name = convertCharsToStrings(session);

        % Plots colrs
        if (strcmp(mouse_info.phenotype,'ctr') || strcmp(mouse_info.phenotype,'WT'))
            clr_m = [74,98,116]./256;
        else
            clr_m = [216,27,96]./256;
        end

        miss_clr = [140,86,86]./256;
        cat1_clr=[36 62 54]./256;
        cat2_clr=[124 169 130]./256;
        cat3_clr=[224 238 198]./256;

        % Single Session
        close all
        disp(strcat('Running session',{' '},sess_name,'...'))
        new_folder_mat = strcat(mat_folder,filesep,group,filesep,setup,filesep,mouse,filesep,sess_name);
        new_folder_out = strcat(out_folder,filesep,group,filesep,setup,filesep,mouse,filesep,char(sess_name),filesep,'reaches_singleLabel_analysis');
        if ~exist("new_folder_out","dir"), mkdir(new_folder_out); end

        %% Variables to plot
        labels={'balistic reaches: ';'reach & hold spout: ';'hidden start: '};

        % IRI
        pks_time_sec = sessions.pks_frames./sessions.video.frame_rate;
        iri_sec = diff(pks_time_sec);

        % Time from ITI to first reach
        % Distributions on time to do first reach
        logic_1st_reach = logical(cat(1,1,diff(reaches.reach_trial)>0));
        time_to_1st_reach_peak = reaches.reach_range_mat(logic_1st_reach,...
            reaches.reach_params.reach_interval.max_reach)./sessions.video.frame_rate;
        time_reach_start = nan(size(reaches.start_forw.start_forw_mat));
        for tt = 1: size(reaches.start_forw.start_forw_mat,1)
            time_reach_start(tt) = reaches.reach_range_mat(tt,reaches.start_forw.start_forw_mat(tt));
        end
        time_to_1st_reach_start = time_reach_start(logic_1st_reach)./sessions.video.frame_rate;

    
        % REACHES
        reaches_all = reaches.reach_mat;
        cat1_reaches = reaches_all(:,:,reaches.cat_reach==1);
        cat2_reaches = reaches_all(:,:,reaches.cat_reach==2);
        cat3_reaches = reaches_all(:,:,reaches.cat_reach==3);
        if isempty(cat2_reaches)
            cat2_reaches = nan(size(cat1_reaches(:,:,1)));
        end

        tm = reaches.time_range;
        reaches_start_stop = reaches_all;
        reaches_start_stop(1:reaches.start_forw.start_forw_mat,:,:) = nan;
        reaches_start_stop(reaches.start_stop.stop_reach_mat:end,:,:) = nan;

        % Speed of REACHES
        speed_reaches_all = cat(1,nan(1,size(reaches_all,2),size(reaches_all,3)),...
            diff(reaches_all));
        speed_cat1_reaches = speed_reaches_all(:,:,reaches.cat_reach==1);
        speed_cat2_reaches = speed_reaches_all(:,:,reaches.cat_reach==2);
        speed_cat3_reaches = speed_reaches_all(:,:,reaches.cat_reach==3);

        % Distribution of durations
        dur_reaches = reaches.start_forw.dur_forw_mat;
        dur_full = reaches.start_stop.dur_start_stop_mat;

        % Distributions oF speed medians
        speed_med = reaches.speed_meds.median_speed_forward_mov_mat;
        speed_med_mid = reaches.speed_meds.median_speed_forward_mov_start_mat;
        speed_max = reaches.speed_meds.max_speed_forward_mov_mat;

        % Distance to water
        water_pos = median(median(trials.water(:,:,:),3,'omitnan'),1,'omitnan');
        reach_to_water = squeeze(sqrt(...
            (reaches_all(:,1,:)-repmat(water_pos(1),[size(reaches_all,1),1,size(reaches_all,3)])).^2 + ...
            (reaches_all(:,2,:)-repmat(water_pos(2),[size(reaches_all,1),1,size(reaches_all,3)])).^2 ...
            ));

        speed_reach_to_water = cat(1,nan(1,size(reach_to_water,2)),diff(reach_to_water));

        % Categories numbers
        perc1 = numel(find(reaches.cat_reach==1));
        perc2 = numel(find(reaches.cat_reach==2));
        perc3 = numel(find(reaches.cat_reach==3));
        perc_cat = [perc1 perc2 perc3];

        % Outcomes numbers
        porpuse_nr = numel(find(reaches.purpose_reach==1));
        noporpuse_nr = numel(find(reaches.purpose_reach==0));
        hit_nr = numel(find(reaches.hit_reach==1));
        miss_nr = numel(find(reaches.hit_reach==0));
        succ_nr =  numel(find(reaches.success_reach==1));
        fail_nr =  numel(find(reaches.success_reach==0));

        % Numbers of categories and outcomes combinations
        nopp_cat1 = numel(intersect(find(reaches.purpose_reach==0),find(reaches.cat_reach==0)));
        nopp_cat2 = numel(intersect(find(reaches.purpose_reach==0),find(reaches.cat_reach==2)));
        nopp_cat3 = numel(intersect(find(reaches.purpose_reach==0),find(reaches.cat_reach==3)));
        pp_cat1 = numel(intersect(find(reaches.purpose_reach==1),find(reaches.cat_reach==1)));
        pp_cat2 = numel(intersect(find(reaches.purpose_reach==1),find(reaches.cat_reach==2)));
        pp_cat3 = numel(intersect(find(reaches.purpose_reach==1),find(reaches.cat_reach==3)));

        %hit_no = numel(find(reaches.hit_reach==0));
        nohit_cat1 = numel(intersect(find(reaches.hit_reach==0),find(reaches.cat_reach==1)));
        nohit_cat2 = numel(intersect(find(reaches.hit_reach==0),find(reaches.cat_reach==2)));
        nohit_cat3 = numel(intersect(find(reaches.hit_reach==0),find(reaches.cat_reach==3)));
        hit_cat1 = numel(intersect(find(reaches.hit_reach==1),find(reaches.cat_reach==1)));
        hit_cat2 = numel(intersect(find(reaches.hit_reach==1),find(reaches.cat_reach==2)));
        hit_cat3 = numel(intersect(find(reaches.hit_reach==1),find(reaches.cat_reach==3)));

        %suc_no = numel(find(reaches.success_reach==0));
        nosuc_cat1 = numel(intersect(find(reaches.success_reach==0),find(reaches.cat_reach==1)));
        nosuc_cat2 = numel(intersect(find(reaches.success_reach==0),find(reaches.cat_reach==2)));
        nosuc_cat3 = numel(intersect(find(reaches.success_reach==10),find(reaches.cat_reach==3)));
        suc_cat1 = numel(intersect(find(reaches.success_reach==1),find(reaches.cat_reach==1)));
        suc_cat2 = numel(intersect(find(reaches.success_reach==1),find(reaches.cat_reach==2)));
        suc_cat3 = numel(intersect(find(reaches.success_reach==1),find(reaches.cat_reach==3)));

        save(strcat(new_folder_mat,filesep,'reaches_analysis.mat'));

        %% FIGURES
        axeOpt = {'linewidth',1.5,'box','off','GridAlpha',0,'ticklength',[1,1]*.01};

        % IRI
        b_w = 0.1;
        figure(1)
        histogram(iri_sec,'binwidth',b_w,'Normalization','pdf','facecolor',clr_m,'FaceAlpha',0.5);
        xlim([-.5 4])
        xlabel('time (s)'); ylabel('pdf')
        title({'inter-reach interval'...
            sprintf('%i%s%i%s',length(pks_time_sec),' reaches in ',trials.n_trials,' trials')});
        set(gca,axeOpt{:});
        if save_flag, saveas(gcf,strcat(new_folder_out,filesep,'iri.png'),'png'); end


        %%
        b_w = 0.1;
        figure()

        subplot(121)
        histogram(time_to_1st_reach_start,'binwidth',b_w,'Normalization','pdf','facecolor',clr_m,'FaceAlpha',0.5);
        set(gca,axeOpt{:})
        xlabel('time to start 1st reach'); ylabel('pdf')
        xlim([-.5 4])

        subplot(122)
        histogram(time_to_1st_reach_peak,'binwidth',b_w,'Normalization','pdf','facecolor',clr_m,'FaceAlpha',0.5);
        set(gca,axeOpt{:})
        xlabel('time to peak 1st reach'); ylabel('pdf')
        xlim([-.5 4])

        set(gcf,'position',[2073         117        1362         487]);
        if save_flag, saveas(gcf,strcat(new_folder_out,filesep,'time_to_1st_reach.png'),'png'); end



        %% Categories splits
        figure()
        explode=[0 1 1];
        pie_cat = pie(perc_cat,explode);

        pText = findobj(pie_cat,'Type','text');
        percentValues = get(pText,'String');
        combinedtxt = strcat(labels,percentValues);
        pText(1).String = combinedtxt(1);
        pText(2).String = combinedtxt(2);
        pText(3).String = combinedtxt(3);
        % Create legend
        legend(labels,'Location','southwest');
        set(gcf, 'Position', [300 230 700 567])
        if save_flag, saveas(gcf,strcat(new_folder_out,filesep,'cat_pie_chart.png'),'png'); end


        %%  Divede purpose, hit and success dived by cat
        %pp_no = numel(find(reaches.purpose_reach==0));
        fHand = figure(3);
        aHand = axes('parent', fHand);
        bar_class_cat = [...
            %    pp_no 0 0; ...
            nopp_cat1 nopp_cat2 nopp_cat3; ...
            pp_cat1 pp_cat2 pp_cat3; 0 0 0; ...
            %    hit_no 0 0; ...
            nohit_cat1 nohit_cat2 nohit_cat3; ...
            hit_cat1 hit_cat2 hit_cat3; 0 0 0; ...
            %    suc_no 0 0; ...
            nosuc_cat1 nosuc_cat2 nosuc_cat3; ...
            suc_cat1 suc_cat2 suc_cat3];
        b1 = bar(bar_class_cat, 'stacked','FaceColor','flat');
        for k = 1:size(bar_class_cat,2)
            b1(k).CData = k;
        end
        hold(aHand, 'on')
        legend(labels,'Location','northeast');
        set(gca, 'XTick', [1,2,4,5,7,8], 'XTickLabel', {'(no)','purpose','(no)','hit','(no)','success'})
        shg
        hold(aHand, 'off')
        if save_flag, saveas(gcf,strcat(new_folder_out,filesep,'purpose_hit_sucess.png'),'png'); end

        %% Reaches over time
        transp = .05;
        lw_th = 1;
        mw_m = 2;
        x_lim = [0 350];
        z_lim = [100 350];

        figure(4)
        subplot(221)
        plot(tm,squeeze(reaches_all(:,1,:)),'linewidth',lw_th,'color',[0.5 0.5 0.5 transp]); hold on
        plot(tm,median(reaches_all(:,1,:),3,'omitnan'),'linewidth',mw_m,'color',[0 0 0]); hold off
        xlim([tm(1),tm(end)]); ylim(x_lim);
        xlabel('time (s)'); ylabel('x (px)')
        set(gca,axeOpt{:})

        subplot(222)
        plot(tm,squeeze(reaches_all(:,2,:)),'linewidth',lw_th,'color',[0.5 0.5 0.5 transp]); hold on
        plot(tm,median(reaches_all(:,2,:),3,'omitnan'),'linewidth',mw_m,'color',[0 0 0]); hold off
        xlim([tm(1) tm(end)]); ylim(z_lim);
        xlabel('time (s)'); ylabel('z (px)')
        set(gca,axeOpt{:})

        subplot(223)
        plot(tm,squeeze(cat3_reaches(:,1,:)),'linewidth',lw_th,'color',cat(2,cat3_clr,transp)); hold on
        plot(tm,squeeze(cat2_reaches(:,1,:)),'linewidth',lw_th,'color',cat(2,cat2_clr,transp)); hold on
        plot(tm,squeeze(cat1_reaches(:,1,:)),'linewidth',lw_th,'color',cat(2,cat1_clr,transp)); hold on
        plot(tm,median(cat3_reaches(:,1,:),3,'omitnan'),'linewidth',mw_m,'color',cat3_clr);
        plot(tm,median(cat2_reaches(:,1,:),3,'omitnan'),'linewidth',mw_m,'color',cat2_clr);
        plot(tm,median(cat1_reaches(:,1,:),3,'omitnan'),'linewidth',mw_m,'color',cat1_clr);
        xlim([tm(1) tm(end)]); ylim(x_lim);
        xlabel('time (s)'); ylabel('x (px)')
        set(gca,axeOpt{:})

        subplot(224)
        plot(tm,squeeze(cat3_reaches(:,2,:)),'linewidth',lw_th,'color',cat(2,cat3_clr,transp)); hold on
        plot(tm,squeeze(cat2_reaches(:,2,:)),'linewidth',lw_th,'color',cat(2,cat2_clr,transp)); hold on
        plot(tm,squeeze(cat1_reaches(:,2,:)),'linewidth',lw_th,'color',cat(2,cat1_clr,transp)); hold on
        plot(tm,median(cat3_reaches(:,2,:),3,'omitnan'),'linewidth',mw_m,'color',cat3_clr);
        plot(tm,median(cat2_reaches(:,2,:),3,'omitnan'),'linewidth',mw_m,'color',cat2_clr);
        plot(tm,median(cat1_reaches(:,2,:),3,'omitnan'),'linewidth',mw_m,'color',cat1_clr);
        xlim([tm(1) tm(end)]); ylim(z_lim);
        xlabel('time (s)'); ylabel('z (px)')
        set(gca,axeOpt{:})

        set(gcf,'color','w','position',[2395 262 1304 653])
        if save_flag, saveas(gcf,strcat(new_folder_out,filesep,'x_y_overtime.png'),'png'); end

        %% Reaches over time just medians
        figure(5)
        subplot(221)
        plot(tm,squeeze(reaches_all(:,1,:)),'linewidth',lw_th,'color',[0.5 0.5 0.5 transp]); hold on
        plot(tm,median(reaches_all(:,1,:),3,'omitnan'),'linewidth',mw_m,'color',[0 0 0]); hold off
        xlim([tm(1),tm(end)]); ylim(x_lim);
        xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
        xlabel('time (s)'); ylabel('x (px)')
        set(gca,axeOpt{:})

        subplot(222)
        plot(tm,squeeze(reaches_all(:,2,:)),'linewidth',lw_th,'color',[0.5 0.5 0.5 transp]); hold on
        plot(tm,median(reaches_all(:,2,:),3,'omitnan'),'linewidth',mw_m,'color',[0 0 0]); hold off
        xlim([tm(1) tm(end)]); ylim(z_lim);
        xlabel('time (s)'); ylabel('z (px)')
        xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
        set(gca,axeOpt{:})

        subplot(223)

        plot(tm,median(cat3_reaches(:,1,:),3,'omitnan'),'linewidth',mw_m,'color',cat3_clr); hold on
        plot(tm,median(cat2_reaches(:,1,:),3,'omitnan'),'linewidth',mw_m,'color',cat2_clr);
        plot(tm,median(cat1_reaches(:,1,:),3,'omitnan'),'linewidth',mw_m,'color',cat1_clr);
        xlim([tm(1) tm(end)]); ylim(x_lim);
        xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
        xlabel('time (s)'); ylabel('x (px)')
        set(gca,axeOpt{:})

        subplot(224)

        plot(tm,median(cat3_reaches(:,2,:),3,'omitnan'),'linewidth',mw_m,'color',cat3_clr); hold on
        plot(tm,median(cat2_reaches(:,2,:),3,'omitnan'),'linewidth',mw_m,'color',cat2_clr);
        plot(tm,median(cat1_reaches(:,2,:),3,'omitnan'),'linewidth',mw_m,'color',cat1_clr);
        xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
        xlim([tm(1) tm(end)]); ylim(z_lim);
        xlabel('time (s)'); ylabel('z (px)')
        set(gca,axeOpt{:})

        set(gcf,'color','w','position',[2395 262 1304 653])
        if save_flag, saveas(gcf,strcat(new_folder_out,filesep,'x_y_overtime_noSingle.png'),'png'); end

        %% Partial speeds
        ylim_speed = [-20 30];

        figure(6)
        subplot(221)
        plot(tm,squeeze(speed_reaches_all(:,1,:)),'linewidth',lw_th,'color',[0.5 0.5 0.5 transp]); hold on
        plot(tm,median(speed_reaches_all(:,1,:),3,'omitnan'),'linewidth',mw_m,'color',[0 0 0]); hold off
        xlim([tm(1),tm(end)]); ylim(ylim_speed);
        xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
        xlabel('time (s)'); ylabel('speed x (px/s)')
        set(gca,axeOpt{:})

        subplot(222)
        plot(tm,squeeze(speed_reaches_all(:,2,:)),'linewidth',lw_th,'color',[0.5 0.5 0.5 transp]); hold on
        plot(tm,median(speed_reaches_all(:,2,:),3,'omitnan'),'linewidth',mw_m,'color',[0 0 0]); hold off
        xlim([tm(1) tm(end)]); ylim(ylim_speed);
        xlabel('time (s)'); ylabel('speed z (px/s)')
        xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
        set(gca,axeOpt{:})

        subplot(223)

        plot(tm,median(speed_cat3_reaches(:,1,:),3,'omitnan'),'linewidth',mw_m,'color',cat3_clr); hold on
        plot(tm,median(speed_cat2_reaches(:,1,:),3,'omitnan'),'linewidth',mw_m,'color',cat2_clr);
        plot(tm,median(speed_cat1_reaches(:,1,:),3,'omitnan'),'linewidth',mw_m,'color',cat1_clr);
        xlim([tm(1) tm(end)]); %ylim(x_lim);
        xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
        xlabel('time (s)'); ylabel('speed x (px/s)')
        set(gca,axeOpt{:})

        subplot(224)

        plot(tm,median(speed_cat3_reaches(:,2,:),3,'omitnan'),'linewidth',mw_m,'color',cat3_clr); hold on
        plot(tm,median(speed_cat2_reaches(:,2,:),3,'omitnan'),'linewidth',mw_m,'color',cat2_clr);
        plot(tm,median(speed_cat1_reaches(:,2,:),3,'omitnan'),'linewidth',mw_m,'color',cat1_clr);
        xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
        xlim([tm(1) tm(end)]); %ylim(z_lim);
        xlabel('time (s)'); ylabel('speed z (px/s)')
        legend(flipud(labels),'box','off')
        set(gca,axeOpt{:})

        set(gcf,'color','w','position',[2395 262 1304 653])
        if save_flag, saveas(gcf,strcat(new_folder_out,filesep,'speed_x_y_overtime.png'),'png'); end



        %% distance to water

        ylim_speed_dist = [-30 20];

        figure(7)
        subplot(221)
        plot(tm,reach_to_water,'linewidth',lw_th,'color',[0.5 0.5 0.5 transp]); hold on
        plot(tm,median(reach_to_water,2,'omitnan'),'linewidth',mw_m,'color',[0 0 0]); hold off
        xlim([tm(1),tm(end)]); %ylim(x_lim);
        xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
        xlabel('time (s)'); ylabel('ditance to water (px)')
        set(gca,axeOpt{:})
        title('distance to water')

        subplot(223)
        plot(tm,speed_reach_to_water,'linewidth',lw_th,'color',[0.5 0.5 0.5 transp]); hold on
        plot(tm,median(speed_reach_to_water,2,'omitnan'),'linewidth',mw_m,'color',[0 0 0]); hold off
        xlim([tm(1),tm(end)]); ylim(ylim_speed_dist);
        xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
        xlabel('time (s)'); ylabel('speed to water (px)')
        set(gca,axeOpt{:})
        title('speed in relation to water')


        subplot(222)
        plot(tm,median(reach_to_water(:,reaches.cat_reach == 3),2,'omitnan'),'linewidth',mw_m,'color',cat3_clr); hold on
        plot(tm,median(reach_to_water(:,reaches.cat_reach == 2),2,'omitnan'),'linewidth',mw_m,'color',cat2_clr); hold on
        plot(tm,median(reach_to_water(:,reaches.cat_reach == 1),2,'omitnan'),'linewidth',mw_m,'color',cat1_clr); hold on
        xlabel('time (s)'); ylabel('ditance to water (px)')
        xlim([tm(1),tm(end)]); %ylim(x_lim);
        xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
        set(gca,axeOpt{:})



        subplot(224)
        plot(tm,median(speed_reach_to_water(:,reaches.cat_reach == 3),2,'omitnan'),'linewidth',mw_m,'color',cat3_clr); hold on
        plot(tm,median(speed_reach_to_water(:,reaches.cat_reach == 2),2,'omitnan'),'linewidth',mw_m,'color',cat2_clr); hold on
        plot(tm,median(speed_reach_to_water(:,reaches.cat_reach == 1),2,'omitnan'),'linewidth',mw_m,'color',cat1_clr); hold on
        xlabel('time (s)'); ylabel('speed to water (px)')
        xlim([tm(1),tm(end)]); ylim(ylim_speed_dist);
        xline(0,'--','LineWidth',2,'color',[.9 .9 .9])
        set(gca,axeOpt{:})
        legend(flipud(labels(1:3)),'box','off')


        set(gcf,'color','w','position',[2176         176        1521         676])
        if save_flag, saveas(gcf,strcat(new_folder_out,filesep,'distance_speed_to_water.png'),'png'); end






        %% Reaches trajectoires
        max_reach = reaches.reach_params.reach_interval.max_reach;
        rising = reaches.start_forw.start_forw_mat;
        rising_med = round(median(rising));

        axis_px = [0 310 80 350];
        range = 1:size(reaches_all,1);
        range_pt = 1:max_reach;

        figure(8)
        subplot(221)
        plot(squeeze(reaches_all(range,1,:)),squeeze(reaches_all(range,2,:)),'linewidth',lw_th,'color',[0.5 0.5 0.5 transp]);hold on
        plot(median(reaches_all(range,1,:),3,'omitnan'),median(reaches_all(range,2,:),3,'omitnan'),'linewidth',mw_m,'color',[0 0 0]);
        plot(median(reaches_all(max_reach,1,:),3,'omitnan'),median(reaches_all(max_reach,2,:),3,'omitnan'),'ko','LineWidth',5);
        plot(median(reaches_all(rising_med,1,:),3,'omitnan'),median(reaches_all(rising_med,2,:),3,'omitnan'),'ko','LineWidth',1);  hold off
        axis([axis_px])
        set(gca,axeOpt{:})
        xlabel('x'); ylabel('z')

        subplot(222)
        plot(median(cat3_reaches(max_reach,1,:),3,'omitnan'),median(cat3_reaches(max_reach,2,:),3,'omitnan'),'o','LineWidth',5,'color',cat3_clr); hold on
        plot(median(cat3_reaches(range,1,:),3,'omitnan'),median(cat3_reaches(range,2,:),3,'omitnan'),'linewidth',mw_m,'color',cat3_clr);
        plot(median(cat2_reaches(max_reach,1,:),3,'omitnan'),median(cat2_reaches(max_reach,2,:),3,'omitnan'),'o','LineWidth',5,'color',cat2_clr);
        plot(median(cat2_reaches(range,1,:),3,'omitnan'),median(cat2_reaches(range,2,:),3,'omitnan'),'linewidth',mw_m,'color',cat2_clr);
        plot(median(cat1_reaches(max_reach,1,:),3,'omitnan'),median(cat1_reaches(max_reach,2,:),3,'omitnan'),'o','LineWidth',5,'color',cat1_clr);
        plot(median(cat1_reaches(range,1,:),3,'omitnan'),median(cat1_reaches(range,2,:),3,'omitnan'),'linewidth',mw_m,'color',cat1_clr);
        axis([axis_px])
        hold off
        set(gca,axeOpt{:})
        xlabel('x'); ylabel('z')

        subplot(223)
        plot(squeeze(reaches_all(range_pt,1,:)),squeeze(reaches_all(range_pt,2,:)),'linewidth',lw_th,'color',[0.5 0.5 0.5 transp]);hold on
        plot(median(reaches_all(range_pt,1,:),3,'omitnan'),median(reaches_all(range_pt,2,:),3,'omitnan'),'linewidth',mw_m,'color',[0 0 0]);
        plot(median(reaches_all(max_reach,1,:),3,'omitnan'),median(reaches_all(max_reach,2,:),3,'omitnan'),'ko','LineWidth',5);
        plot(median(reaches_all(rising_med,1,:),3,'omitnan'),median(reaches_all(rising_med,2,:),3,'omitnan'),'ko','LineWidth',1);  hold off
        axis([axis_px])
        set(gca,axeOpt{:})
        xlabel('x'); ylabel('z')

        subplot(224)
        plot(median(cat3_reaches(range_pt,1,:),3,'omitnan'),median(cat3_reaches(range_pt,2,:),3,'omitnan'),'linewidth',mw_m,'color',cat3_clr);  hold on
        plot(median(cat2_reaches(range_pt,1,:),3,'omitnan'),median(cat2_reaches(range_pt,2,:),3,'omitnan'),'linewidth',mw_m,'color',cat2_clr);
        plot(median(cat1_reaches(range_pt,1,:),3,'omitnan'),median(cat1_reaches(range_pt,2,:),3,'omitnan'),'linewidth',mw_m,'color',cat1_clr);
        plot(median(cat3_reaches(max_reach,1,:),3,'omitnan'),median(cat3_reaches(max_reach,2,:),3,'omitnan'),'o','LineWidth',5,'color',cat3_clr);
        plot(median(cat2_reaches(max_reach,1,:),3,'omitnan'),median(cat2_reaches(max_reach,2,:),3,'omitnan'),'o','LineWidth',5,'color',cat2_clr);
        plot(median(cat1_reaches(max_reach,1,:),3,'omitnan'),median(cat1_reaches(max_reach,2,:),3,'omitnan'),'o','LineWidth',5,'color',cat1_clr);
        axis([axis_px])
        legend(flipud(labels),'Location','southeast','box', 'off');
        hold off
        set(gca,axeOpt{:})
        xlabel('x'); ylabel('z')

        set(gcf,'color','w','position',[2405         138        963         825]);
        if save_flag, saveas(gcf,strcat(new_folder_out,filesep,'reaches_xz.png'),'png'); end

        %%
        nbins = 50;
        figure()
        subplot(221)
        plot(squeeze(reaches_start_stop(:,1,:)),squeeze(reaches_start_stop(:,2,:)),'linewidth',lw_th,'color',[0.5 0.5 0.5 transp]);hold on
        plot(median(reaches_start_stop(:,1,:),3,'omitnan'),median(reaches_start_stop(:,2,:),3,'omitnan'),'linewidth',mw_m,'color',[0 0 0]);
        plot(median(reaches_start_stop(max_reach,1,:),3,'omitnan'),median(reaches_start_stop(max_reach,2,:),3,'omitnan'),'ko','LineWidth',5);
        axis([axis_px])
        set(gca,axeOpt{:})
        xlabel('x'); ylabel('z')
        title('reaches start stop')

        subplot(222)
        plot(median(reaches_start_stop(:,1,reaches.cat_reach==3),3,'omitnan'),median(reaches_start_stop(:,2,reaches.cat_reach==3),3,'omitnan'),'linewidth',mw_m,'color',cat3_clr); hold on
        plot(median(reaches_start_stop(:,1,reaches.cat_reach==1),3,'omitnan'),median(reaches_start_stop(:,2,reaches.cat_reach==1),3,'omitnan'),'linewidth',mw_m,'color',cat1_clr);
        plot(median(reaches_start_stop(:,1,reaches.cat_reach==2),3,'omitnan'),median(reaches_start_stop(:,2,reaches.cat_reach==2),3,'omitnan'),'linewidth',mw_m,'color',cat2_clr);
        plot(median(reaches_start_stop(max_reach,1,reaches.cat_reach==3),3,'omitnan'),median(reaches_start_stop(max_reach,2,reaches.cat_reach==3),3,'omitnan'),'o','LineWidth',5,'color',cat3_clr);
        plot(median(reaches_start_stop(max_reach,1,reaches.cat_reach==2),3,'omitnan'),median(reaches_start_stop(max_reach,2,reaches.cat_reach==2),3,'omitnan'),'o','LineWidth',5,'color',cat2_clr);
        plot(median(reaches_start_stop(max_reach,1,reaches.cat_reach==1),3,'omitnan'),median(reaches_start_stop(max_reach,2,reaches.cat_reach==1),3,'omitnan'),'o','LineWidth',5,'color',cat1_clr);
        axis([axis_px])
        set(gca,axeOpt{:})
        xlabel('x'); ylabel('z')
        legend(flipud(labels),'Location','southeast','box', 'off');
        title('reaches start stop categories')

        subplot(223)
        histogram2(squeeze(reaches_all(:,1,:)),...
            squeeze(reaches_all(:,2,:)),...
            nbins,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none',...
            'normalization','pdf');
        set(gca,'CLim',[0 7E-5])
        colorbar
        axis([axis_px])
        set(gca,axeOpt{:})
        xlabel('x'); ylabel('z')
        title('full reach')

        subplot(224)
        histogram2(squeeze(reaches_start_stop(:,1,:)),...
            squeeze(reaches_start_stop(:,2,:)),...
            nbins,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none',...
            'normalization','pdf');
        set(gca,'CLim',[0 7E-5])
        colorbar
        axis([axis_px])
        set(gca,axeOpt{:})
        xlabel('x'); ylabel('z')
        title('reach from start to stop')

        set(gcf,'color','w','position',[2405         138        963         825]);
        if save_flag, saveas(gcf,strcat(new_folder_out,filesep,'reaches_dist.png'),'png'); end


        %% Reaches categories by outcome
        figure(9)
        range = range_pt;

        subplot(131)
        trials_class_y = reaches.purpose_reach==1;
        trials_class_n = reaches.purpose_reach==0;

        plot(squeeze(reaches_all(range,1,trials_class_y)),squeeze(reaches_all(range,2,trials_class_y)),'linewidth',lw_th,'color',[0.5 0.5 0.5 transp]);hold on
        plot(squeeze(reaches_all(range,1,trials_class_n)),squeeze(reaches_all(range,2,trials_class_n)),'linewidth',lw_th,'color',cat(2,miss_clr,transp));
        plot(median(reaches_all(range,1,trials_class_y),3,'omitnan'),median(reaches_all(range,2,trials_class_y),3,'omitnan'),'linewidth',mw_m,'color',[0 0 0]);
        plot(median(reaches_all(range,1,trials_class_n),3,'omitnan'),median(reaches_all(range,2,trials_class_n),3,'omitnan'),'linewidth',mw_m,'color',miss_clr);
        plot(median(reaches_all(max_reach,1,trials_class_y),3,'omitnan'),median(reaches_all(max_reach,2,trials_class_y),3,'omitnan'),'ko','LineWidth',5);
        plot(median(reaches_all(max_reach,1,trials_class_n),3,'omitnan'),median(reaches_all(max_reach,2,trials_class_n),3,'omitnan'),'o','LineWidth',5,'color',miss_clr);
        axis([axis_px])
        set(gca,axeOpt{:})
        xlabel('x'); ylabel('z')
        title('purposefull vs not')

        subplot(132)
        trials_class_y = reaches.hit_reach==1;
        trials_class_n = reaches.hit_reach==0;

        plot(squeeze(reaches_all(range,1,trials_class_y)),squeeze(reaches_all(range,2,trials_class_y)),'linewidth',lw_th,'color',[0.5 0.5 0.5 transp]);hold on
        plot(squeeze(reaches_all(range,1,trials_class_n)),squeeze(reaches_all(range,2,trials_class_n)),'linewidth',lw_th,'color',cat(2,miss_clr,transp));
        plot(median(reaches_all(range,1,trials_class_y),3,'omitnan'),median(reaches_all(range,2,trials_class_y),3,'omitnan'),'linewidth',mw_m,'color',[0 0 0]);
        plot(median(reaches_all(range,1,trials_class_n),3,'omitnan'),median(reaches_all(range,2,trials_class_n),3,'omitnan'),'linewidth',mw_m,'color',miss_clr);
        plot(median(reaches_all(max_reach,1,trials_class_y),3,'omitnan'),median(reaches_all(max_reach,2,trials_class_y),3,'omitnan'),'ko','LineWidth',5);
        plot(median(reaches_all(max_reach,1,trials_class_n),3,'omitnan'),median(reaches_all(max_reach,2,trials_class_n),3,'omitnan'),'o','LineWidth',5,'color',miss_clr);
        axis([axis_px])
        set(gca,axeOpt{:})
        xlabel('x'); ylabel('z')
        title('hit vs miss')

        subplot(133)
        trials_class_y = reaches.success_reach==1;
        trials_class_n = reaches.success_reach==0;

        plot(squeeze(reaches_all(range,1,trials_class_y)),squeeze(reaches_all(range,2,trials_class_y)),'linewidth',lw_th,'color',[0.5 0.5 0.5 transp]);hold on
        plot(squeeze(reaches_all(range,1,trials_class_n)),squeeze(reaches_all(range,2,trials_class_n)),'linewidth',lw_th,'color',cat(2,miss_clr,transp));
        plot(median(reaches_all(range,1,trials_class_y),3,'omitnan'),median(reaches_all(range,2,trials_class_y),3,'omitnan'),'linewidth',mw_m,'color',[0 0 0]);
        plot(median(reaches_all(range,1,trials_class_n),3,'omitnan'),median(reaches_all(range,2,trials_class_n),3,'omitnan'),'linewidth',mw_m,'color',miss_clr);
        plot(median(reaches_all(max_reach,1,trials_class_y),3,'omitnan'),median(reaches_all(max_reach,2,trials_class_y),3,'omitnan'),'ko','LineWidth',5);
        plot(median(reaches_all(max_reach,1,trials_class_n),3,'omitnan'),median(reaches_all(max_reach,2,trials_class_n),3,'omitnan'),'o','LineWidth',5,'color',miss_clr);
        axis([axis_px])
        set(gca,axeOpt{:})
        xlabel('x'); ylabel('z')
        title('success vs faill')

        set(gcf,'color','w','position',[ 1996         298        1839         534]);
        if save_flag, saveas(gcf,strcat(new_folder_out,filesep,'outcome_reaches.png'),'png'); end


        %% Distribtuion of speed from start to max reach
        b_w_s = 2;
        axis_speed_h = [-50 50];

        figure(10)
        subplot(221)
        histogram(speed_reaches_all(1:max_reach,1,:),'binwidth',b_w_s,'Normalization','pdf','facecolor','k','FaceAlpha',0.5); hold on
        xlim(axis_speed_h);
        set(gca,axeOpt{:})
        xlabel('speed reach forward in x (px/s)'); ylabel('counts')

        subplot(222)
        histogram(speed_reaches_all(1:max_reach,2,:),'binwidth',b_w_s,'Normalization','pdf','facecolor','k','FaceAlpha',0.5); hold on
        xlim(axis_speed_h);
        set(gca,axeOpt{:})
        xlabel('speed reach forward in z (px/s)'); ylabel('counts')

        subplot(223)
        histogram(speed_reaches_all(1:max_reach,1,reaches.cat_reach==3),'binwidth',b_w_s,'Normalization','pdf','facecolor',cat3_clr,'FaceAlpha',0.5); hold on
        histogram(speed_reaches_all(1:max_reach,1,reaches.cat_reach==2),'binwidth',b_w_s,'Normalization','pdf','facecolor',cat2_clr,'FaceAlpha',0.5); hold on
        histogram(speed_reaches_all(1:max_reach,1,reaches.cat_reach==1),'binwidth',b_w_s,'Normalization','pdf','facecolor',cat1_clr,'FaceAlpha',0.5); hold on
        xlim(axis_speed_h);
        set(gca,axeOpt{:})
        xlabel('speed reach forward in x (px/s)'); ylabel('counts')

        subplot(224)
        histogram(speed_reaches_all(1:max_reach,2,reaches.cat_reach==3),'binwidth',b_w_s,'Normalization','pdf','facecolor',cat3_clr,'FaceAlpha',0.5); hold on
        histogram(speed_reaches_all(1:max_reach,2,reaches.cat_reach==2),'binwidth',b_w_s,'Normalization','pdf','facecolor',cat2_clr,'FaceAlpha',0.5); hold on
        histogram(speed_reaches_all(1:max_reach,2,reaches.cat_reach==1),'binwidth',b_w_s,'Normalization','pdf','facecolor',cat1_clr,'FaceAlpha',0.5); hold on
        xlim(axis_speed_h);
        set(gca,axeOpt{:})
        xlabel('speed reach forward in z (px/s)'); ylabel('counts')
        set(gcf,'color','w','position',[2346          84        1137         892])
        if save_flag, saveas(gcf,strcat(new_folder_out,filesep,'speed_forward_distributions.png'),'png'); end

        %% Distribution of durations

        figure(11)
        subplot(221)
        histogram(dur_reaches,'binwidth',0.01,'Normalization','pdf','facecolor','k','FaceAlpha',0.5); hold on
        set(gca,axeOpt{:})
        xlabel('duration forward (s)'); ylabel('counts')

        subplot(223)
        histogram(dur_reaches(reaches.cat_reach==2),'binwidth',0.01,'Normalization','pdf','facecolor',cat2_clr,'FaceAlpha',0.5); hold on
        histogram(dur_reaches(reaches.cat_reach==1),'binwidth',0.01,'Normalization','pdf','facecolor',cat1_clr,'FaceAlpha',0.5); hold on
        set(gca,axeOpt{:})
        xlabel('duration forward (s)'); ylabel('counts')

        subplot(222)
        histogram(dur_full,'binwidth',0.01,'Normalization','pdf','facecolor','k','FaceAlpha',0.5); hold on
        set(gca,axeOpt{:})
        xlabel('duration start-stop (s)'); ylabel('counts')

        subplot(224)
        histogram(dur_full(reaches.cat_reach==2),'binwidth',0.01,'Normalization','pdf','facecolor',cat2_clr,'FaceAlpha',0.5); hold on
        histogram(dur_full(reaches.cat_reach==1),'binwidth',0.01,'Normalization','pdf','facecolor',cat1_clr,'FaceAlpha',0.5); hold on
        set(gca,axeOpt{:})
        xlabel('duration start-stop (s)'); ylabel('counts')
        legend(flipud(labels(1:2)),'box','off')
        set(gcf,'color','w','position',[2346          84        1137         892])
        if save_flag, saveas(gcf,strcat(new_folder_out,filesep,'durations.png'),'png'); end


        %% Distributions on speed medians
        b_w = 1;
        figure(12)

        subplot(231)
        histogram(speed_med,'binwidth',b_w,'Normalization','pdf','facecolor','k','FaceAlpha',0.5); hold on
        set(gca,axeOpt{:})
        xlabel('median speed (px/s)'); ylabel('counts')
        xlim([-5 25])

        subplot(234)
        histogram(speed_med(reaches.cat_reach==3),'binwidth',b_w,'Normalization','pdf','facecolor',cat3_clr,'FaceAlpha',0.5); hold on
        histogram(speed_med(reaches.cat_reach==2),'binwidth',b_w,'Normalization','pdf','facecolor',cat2_clr,'FaceAlpha',0.5); hold on
        histogram(speed_med(reaches.cat_reach==1),'binwidth',b_w,'Normalization','pdf','facecolor',cat1_clr,'FaceAlpha',0.5); hold on
        set(gca,axeOpt{:})
        xlim([-5 25])
        xlabel('median speed  (px/s)'); ylabel('counts')

        subplot(232)
        histogram(speed_med_mid,'binwidth',b_w,'Normalization','pdf','facecolor','k','FaceAlpha',0.5); hold on
        set(gca,axeOpt{:})
        xlim([-5 25])
        xlabel('median speed midreach (px/s)'); ylabel('counts')

        subplot(235)
        histogram(speed_med_mid(reaches.cat_reach==3),'binwidth',b_w,'Normalization','pdf','facecolor',cat3_clr,'FaceAlpha',0.5); hold on
        histogram(speed_med_mid(reaches.cat_reach==2),'binwidth',b_w,'Normalization','pdf','facecolor',cat2_clr,'FaceAlpha',0.5); hold on
        histogram(speed_med_mid(reaches.cat_reach==1),'binwidth',b_w,'Normalization','pdf','facecolor',cat1_clr,'FaceAlpha',0.5); hold on
        set(gca,axeOpt{:})
        xlabel('median speed midreach (px/s)'); ylabel('counts')
        xlim([-5 25])

        subplot(233)
        histogram(speed_max,'binwidth',b_w,'Normalization','pdf','facecolor','k','FaceAlpha',0.5); hold on
        set(gca,axeOpt{:})
        xlabel('max speed (px/ss)'); ylabel('counts')
        xlim([-5 50])

        subplot(236)
        histogram(speed_max(reaches.cat_reach==3),'binwidth',b_w,'Normalization','pdf','facecolor',cat3_clr,'FaceAlpha',0.5); hold on
        histogram(speed_max(reaches.cat_reach==2),'binwidth',b_w,'Normalization','pdf','facecolor',cat2_clr,'FaceAlpha',0.5); hold on
        histogram(speed_max(reaches.cat_reach==1),'binwidth',b_w,'Normalization','pdf','facecolor',cat1_clr,'FaceAlpha',0.5); hold on
        set(gca,axeOpt{:})
        xlabel('max speed (px/ss)'); ylabel('counts')
        xlim([-5 50])



        legend(flipud(labels(1:3)),'box','off')
        set(gcf,'color','w','position',[2279         119        1480         857])
        if save_flag, saveas(gcf,strcat(new_folder_out,filesep,'speed_meds.png'),'png'); end



    end
end


%legend(flipud(labels(1:3)),'box','off')
%set(gcf,'color','w','position',[2279         119        1480         857])
%if save_flag, saveas(gcf,strcat(new_folder_out,filesep,'speed_meds.png'),'png'); end




