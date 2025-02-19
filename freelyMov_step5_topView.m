clear; close all; clc

%% % Paths for each mouse and initializaition variables
raw_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\raw_data';
mat_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\analyzed_data\mat_files';
out_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\analyzed_data\output_files';

%data_pathA = '/Users/teresaserradasduarte/Dropbox (Learning Lab)/Learning Lab Team Folder/Patlab protocols/data/TD/raw_data';
project_name = '20240729_A2aCaspGroup';
setup = 'freely_mov';
datamat_path = strcat(mat_folder,filesep,project_name,filesep,setup);
dataraw_path = strcat(raw_folder,filesep,project_name,filesep,setup);
save_path = strcat(out_folder,filesep,project_name,filesep,setup,filesep,'topView');
if ~exist('save_path','dir'), mkdir(save_path); end

%mouse = 'WT2';
%sess = 'S7';
%raw_path = strcat(raw_folder,filesep,project_name,filesep,setup,filesep,mouse,filesep,sess,filesep);

%% Load
fps=120;
max_time = fps*60*40;
resnet = 'DLC_resnet50_FreelyMovReaching_TopViewAug10shuffle1_950000';
likelihood_threshold = 0.9;

% Find mouse folders
folders_mice = dir(dataraw_path);
folders_mice=folders_mice(3:end,:);
num_animals=size(folders_mice,1);
% Allocate space for mouse cells
mice = cell(num_animals,1);
mice_path = cell(num_animals,1);
num_sessions =3;
% Allocate space for data variables
top_head_all = nan(max_time,2,num_sessions,num_animals);

for m = 1: num_animals
    mice(m,1)=cellstr(convertCharsToStrings(folders_mice(m,1).name));
    mice_path(m,1) = strcat(dataraw_path, filesep, mice(m,1));
    % Find folders (1 per session)
    folders_session = dir(char(mice_path(m,1)));
    folders_session = folders_session(end-2:end,:); % RANGE ADAPTD TO SESISONS OF INTEREST
    num_sessions = size(folders_session,1);
    % Allocate space for session cells
    sessions = cell(num_sessions,1);
    sessions_path = cell(num_sessions,1);

    for s=1:num_sessions
        sessions(s,1)=cellstr(convertCharsToStrings(folders_session(s,1).name));
        sessions_path(s,1) = strcat(dataraw_path, filesep, mice(m,1), filesep, sessions(s,1), filesep);

        % Find csvs with tracking
        searchstr=strcat(resnet,'.csv');
        csvs = wildcardsearch(sessions_path(s,1), searchstr, true);

        if ~isempty(csvs)
            topView=csvread(char(csvs),3,1);

            lik = repmat(topView(:,3),[1, 2]);
            top_head = topView(:,1:2);
            top_head(lik<likelihood_threshold) = nan;

            top_head_all(1:size(top_head,1),:,s,m) = top_head;
        end
    end
end

%%
dt = 1/fps;
speed_all_xy = diff(top_head_all);
speed_all_xy_dt = speed_all_xy./dt;
% speed (velocity vector magnitude)
speed_all_dt = squeeze(sqrt(...
    speed_all_xy_dt(:,1,:,:).^2 + ...
    speed_all_xy_dt(:,2,:,:).^2 ...
    ));


%% Plot stuff
a2a_range = 1:2;
wts_range = 3:4;

a2a_clr = [216,27,96]./256;
wts_clr = [74,98,116]./256;

cat1_clr=[36 62 54]./256;
cat2_clr=[124 169 130]./256;
cat3_clr=[224 238 198]./256;

clrs_m = cat(1,a2a_clr,a2a_clr,wts_clr,wts_clr);


axeOpt = {'linewidth',1.5,'box','off','GridAlpha',0,'ticklength',[1,1]*.01};
axeOpt_H = {'linewidth',1.5,'XTick', [], 'YTick', [],'box','off','GridAlpha',0,};

%%

xy_lims = [200 1000 100 900];
transpa = 0.1;
for sess = 1:num_sessions
    figure
    figg=tiledlayout(2,2);

    for m = 1:num_animals
        nexttile
        scatter(top_head_all(:,1,sess,m),top_head_all(:,2,sess,m),10,'filled',"MarkerFaceColor",clrs_m(m,:),'MarkerFaceAlpha',transpa)
        axis(xy_lims)
        set(gca,axeOpt{:})
        xlabel('x (px)'); ylabel('y (px)')
    end

    title(figg,{'Position on box';char(sessions(sess))})
    set(gcf,'color','w','position',[ 2228         134         937         828])
    saveas(gcf,strcat(save_path,filesep,'position_xy',char(sessions(sess)),'.png'),'png')

end

%%

n_bins = 30;
hist2d_all = nan(n_bins,n_bins,num_sessions,num_animals);
hist2d_all_norm = nan(n_bins,n_bins,num_sessions,num_animals);
entropy = nan(num_sessions,num_animals);

clip_max = 500;
for sess = 1:num_sessions
    figure
    figg=tiledlayout(2,2);

    for m=1:num_animals
        nexttile
        b = histogram2(top_head_all(:,1,sess,m),top_head_all(:,2,sess,m),n_bins,'DisplayStyle','tile','ShowEmptyBins','on','EdgeColor','none');
        axis(xy_lims);
        set(gca,axeOpt_H{:},'CLim',[0 clip_max])
        xlabel('x (px)'); ylabel('y (px)')

        hist2d_all(:,:,sess,m) = b.Values;
        hist2d_all_norm(:,:,sess,m) = b.Values./sum(b.Values(:));

        p = b.Values(:)./sum(b.Values(:));
        p(p==0)=nan;
        entropy(sess,m) = -sum(p.*log2(p),'omitnan');
        nans_tot = length(find(isnan(p)));
        sat_tot = length(find(b.Values(:)>clip_max));

        title({char(mice(m)); sprintf('%s%.3f%s%i',...
            'Entropy = ',entropy(sess,m),' , nr nans = ',nans_tot);...
            sprintf('%s%i','Nr saturated  = ',sat_tot)});

    end
    title(figg,{'Occupancy on box';char(sessions(sess))})
    set(gcf,'color','w','position',[2228          88         937         874])
    saveas(gcf,strcat(save_path,filesep,'occupancy_xy',char(sessions(sess)),'.png'),'png')

end

%%
hist2d_all_nan = hist2d_all;
hist2d_all_nan(hist2d_all==0)=nan;
bw = 5;
for sess = 1:num_sessions
    figure
    figg=tiledlayout(2,2);

    for m=1:num_animals
        nexttile
        histogram(hist2d_all_nan(:,:,sess,m),'BinWidth',bw,'Normalization','pdf','facecolor',clrs_m(m,:))
        %axis(xy_lims);
        set(gca,axeOpt{:})
        xlim([-10 200])
        xlabel('density of occupation'); ylabel('counts')
        title(char(mice(m)),'Color',clrs_m(m,:))
    end
    
    title(figg,{'Distribution of occupancy';char(sessions(sess))})
    set(gcf,'color','w','position',[ 2228         134         937         828])
    saveas(gcf,strcat(save_path,filesep,'occupancy_hist_nan',char(sessions(sess)),'.png'),'png')
end

%%
bw = 5;

hist2d_all_wts = hist2d_all_nan(:,:,2:3,wts_range);
hist2d_all_a2a = hist2d_all_nan(:,:,2:3,a2a_range);
figure
histogram(hist2d_all_wts(~isnan(hist2d_all_wts)),'BinWidth',bw,'Normalization','pdf','facecolor',wts_clr); hold on
histogram(hist2d_all_a2a(~isnan(hist2d_all_a2a)),'BinWidth',bw,'Normalization','pdf','facecolor',a2a_clr); hold on
xlim([-10 200])
set(gca,axeOpt{:})
title('Occupancy density distribution pooled')
xlabel('occupancy density'); ylabel('counts')
saveas(gcf,strcat(save_path,filesep,'speed_hist_occupancy_S6-7.png'),'png')



%%
bw = 30;
for sess = 1:num_sessions
    figure
    figg=tiledlayout(2,2);

    for m=1:num_animals
        nexttile
        histogram(speed_all_dt(:,sess,m),'BinWidth',bw,'Normalization','pdf','facecolor',clrs_m(m,:))
        %axis(xy_lims);
        set(gca,axeOpt{:})
        %xlim([0 1500])
        xlim([1500 3000])
ylim([0 2E-5])
        xlabel('speed (px/s)'); ylabel('counts')
        title(char(mice(m)),'Color',clrs_m(m,:))
    end
    
    title(figg,{'Distribution of speed';char(sessions(sess))})
    set(gcf,'color','w','position',[ 2228         134         937         828])
    saveas(gcf,strcat(save_path,filesep,'speed_hist-zomm',char(sessions(sess)),'.png'),'png')
end

%%
bw = 10;
speed_all_wts = speed_all_dt(:,2:3,wts_range);
speed_all_a2a = speed_all_dt(:,2:3,a2a_range);

figure
subplot(121)
histogram(speed_all_a2a(~isnan(speed_all_a2a)),'BinWidth',bw,'Normalization','pdf','facecolor',a2a_clr); hold on
histogram(speed_all_wts(~isnan(speed_all_wts)),'BinWidth',bw,'Normalization','pdf','facecolor',wts_clr); hold on
xlim([-10 500])
set(gca,axeOpt{:})
xlabel('speed (px/s)'); ylabel('counts')

subplot(122)
histogram(speed_all_a2a(~isnan(speed_all_a2a)),'BinWidth',bw,'Normalization','pdf','facecolor',a2a_clr); hold on
histogram(speed_all_wts(~isnan(speed_all_wts)),'BinWidth',bw,'Normalization','pdf','facecolor',wts_clr); hold on
xlim([1500 2500])
ylim([0 2E-5])
set(gca,axeOpt{:})
title('Speed Distribution pooled')
xlabel('speed (px/s)'); ylabel('counts')
    set(gcf,'color','w','position',[2350         221        1424         518])
saveas(gcf,strcat(save_path,filesep,'speed_hist_pooled_zoom_S6-7.png'),'png')



%%
clip_speed = 1500;
figure();
for sess = 1:num_sessions
    figure
    figg=tiledlayout(2,2);

    for m=1:num_animals
        nexttile
        c = [nan;speed_all_dt(:,sess,m)];
        scatter(top_head_all(:,1,sess,m),top_head_all(:,2,sess,m),20,c,'filled','MarkerFaceAlpha',0.6)
        set(gca, 'XTick', [], 'YTick', [],'Box','on','CLim',[0 clip_speed],axeOpt_H{:})
        axis square;
        xlabel('x (px)'); ylabel('y (px)')
        axis(xy_lims)
        cb = colorbar;
        ylabel(cb,'speed (px/s)')
        title(char(mice(m)),'Color',clrs_m(m,:))


sat_speed = length(find(c>clip_speed));
           

        title({char(mice(m)); ...
            sprintf('%s%i','Nr saturated  = ',sat_speed)});


    end

   title(figg,{'Instantaneous speed at each position';char(sessions(sess))})
    set(gcf,'color','w','position',[ 2228         134         937         828])
    saveas(gcf,strcat(save_path,filesep,'speed_pos',char(sessions(sess)),'.png'),'png')

end











