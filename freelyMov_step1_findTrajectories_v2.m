%% Reaching Freely Moving
clear; close all; clc

%% Load side view
raw_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\raw_data';
mat_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\analyzed_data\mat_files';
out_folder = 'D:\Learning Lab Dropbox\Learning Lab Team Folder\Patlab protocols\Data\TD\behavior_data\analyzed_data\output_files';

group = '20250106_A2aCasp_G2';
setup = 'freely_mov';
mouse = '02_Str34';
%mouse = '02_Str34';
session = 'S2';



all_sess_path_mat=strcat(mat_folder,filesep,group,filesep,setup,filesep,mouse);
if ~exist(all_sess_path_mat,'dir'), mkdir(all_sess_path_mat); end

% Mouse paw preference
% R -> righties | L-> Lefties
% phenotype: ctr or dtx
phenotype = 'a2a';
%paw_pref = 'L';

% Save mouse info
% Mouse info
mouse_info.group = group;
mouse_info.mouse = mouse;
mouse_info.phenotype = phenotype;

% Video / tracking parameters
likelihood_threshold = 0.9;
frame_rate = 120;    % Sampling frequency
T = 1/frame_rate;    % Sampling period
front_tracked_flag = true;

% Display and save
show_plot = false;
show_overview_reach_plot = false;
save_fig_flag = false;

% Figures handles
figOpt = {'color','w'};
axeOpt = {'linewidth',1.5,'box','off','GridAlpha',0.05,'ticklength',[1,1]*.01};

% Find folders (1 per session)
% all_sess_path=strcat(raw_folder,filesep,group,filesep,setup,filesep,mouse);
% folders_session = dir(all_sess_path);
% folders_session = folders_session(3:end,:);
% num_sessions = size(folders_session,1);
% sess = cell(num_sessions,1);
% folders.all_sess_folders = folders_session.name;
% folders.num_sessions = num_sessions;
%save(strcat(all_sess_path_mat,filesep,'mouse_info.mat'),'mouse_info','num_sessions','folders_session','px2mm','')

%% LOOP SESSIONS

%disp('Reading csvs and finding reaches for each session...')

raw_side = strcat(raw_folder, filesep, group, filesep, setup, filesep, mouse, filesep,  session, filesep, 'videosCutSide',filesep,'side');
raw_front = strcat(raw_folder, filesep, group, filesep, setup, filesep, mouse, filesep,  session, filesep, 'videosCutSide',filesep,'front');
out_mat_folder =  strcat(mat_folder, filesep, group, filesep, setup, filesep, mouse, filesep, session,filesep);
output_folder = strcat(out_folder, filesep, group, filesep, setup, filesep, mouse, filesep, session,filesep,'raw_reaches');
if save_fig_flag, if ~exist("output_folder","dir"), mkdir(output_folder); end, end
if ~exist("out_mat_folder","dir"), mkdir(out_mat_folder); end

if  exist(strcat(out_mat_folder,filesep,'session_raw_data.mat'),'file')~=2
    % video params
    VideoMeta = VideoReader(strcat(raw_side,filesep,session,'-camA-trial0.avi'));
    width = VideoMeta.width;
    height = VideoMeta.height;

    VideoMetaFront = VideoReader(strcat(raw_front,filesep,session,'-camB-trial0.avi'));
    widthB = VideoMetaFront.width;
    heightB = VideoMetaFront.height;

    % Data information
    % DNN
    resnet = 'DLC_resnet50_FreelyMovReaching_Side_v2Jan17shuffle1_1000000';
    resnet_front = 'DLC_resnet50_FreelyMovReaching_Front_v2Feb6shuffle1_800000';

    % Find csvs with tracking
    searchstr=strcat(resnet,'.csv');
    csvs = wildcardsearch(raw_side, searchstr, true);

    % Check nr of trials
    indA = find(contains(csvs,'camA'))';
    n_trials = length(indA)-1;
    siz = zeros(n_trials,1);

    % Find max number of frames in 1 trial
    tic
    disp('Loop through trials to find the n_max_points...')
    for i=1:n_trials % to ensure cronologinal ordering
        % cam A
        tmp_camA_csv = strcat(raw_side, filesep,session,'-camA-trial',num2str(i),resnet,'.csv');
        tmp_camA=csvread(tmp_camA_csv,3,1);
        siz(i)=size(tmp_camA,1);
    end
    [max_siz, max_siz_loc]=max(siz); % Check max size to define the mat length
    disp(strcat('Done: max= ',num2str(max_siz)))

    %
    % Allocate space
    n_features_side = 4;
    n_features_front = 4;
    n_max_points = max_siz;
    camA = nan(n_max_points,n_features_side*3,n_trials);
    camB = nan(n_max_points,n_features_front*3,n_trials);

    % time
    time = nan(n_max_points,n_trials);

    % Read from csv
    tic
    disp('finding csvs ordered by cam/view...')
    camA_in_row = [];
    camB_in_row = [];
    for i=1:n_trials % to ensure cronologinal ordering
        % cam A
        tmp_camA_csv = strcat(raw_side, filesep,session,'-camA-trial',num2str(i),resnet,'.csv');
        tmp_camA=csvread(tmp_camA_csv,3,1);
        camA(1:size(tmp_camA,1),:,i) = tmp_camA;
        time(1:size(tmp_camA,1),i) = (1:size(tmp_camA,1))./frame_rate;
        camA_in_row=cat(1,camA_in_row,tmp_camA);

        % cam B
        if front_tracked_flag==true
            tmp_camB_csv = strcat(raw_front, filesep,session,'-camB-trial',num2str(i),resnet_front,'.csv');
            tmp_camB=csvread(tmp_camB_csv,3,1);
            camB(1:size(tmp_camB,1),:,i) = tmp_camB;
            camB_in_row=cat(1,camB_in_row,tmp_camB);
        end
    end
    disp('Done!')
    toc

    % Arrange data - organize in x,y,z by feature (paws and water_A)
    % RIGHT PAW
    % x, z
    tpaw_A = camA(:,1:2,:);
    lik_paw_A = repmat(camA(:,3,:),[1, 2, 1]);
    paw_A_flip_h = abs(tpaw_A(:,2,:)-height-1);
    paw_A=cat(2,tpaw_A(:,1,:),paw_A_flip_h);
    paw_A(lik_paw_A<likelihood_threshold) = nan;
    pawA_in_row = cat(2,camA_in_row(:,1),abs(camA_in_row(:,2)-height-1));

    % SNOUT
    tsnout_A = camA(:,4:5,:);
    lik_snout_A = repmat(camA(:,6,:),[1, 2, 1]);
    snout_A_flip_h = abs(tsnout_A(:,2,:)-height-1);
    snout_A=cat(2,tsnout_A(:,1,:),snout_A_flip_h);
    snout_A(lik_snout_A<likelihood_threshold) = nan;
    snoutA_in_row = cat(2,camA_in_row(:,4),abs(camA_in_row(:,5)-height-1));

    % WATER
    twater_A = camA(:,7:8,:);
    lik_water_A = repmat(camA(:,9,:),[1, 2, 1]);
    water_A_flip_h = abs(twater_A(:,2,:)-height-1);
    water_A=cat(2,twater_A(:,1,:),water_A_flip_h);
    water_A(lik_water_A<likelihood_threshold) = nan;
    waterA_in_row = cat(2,camA_in_row(:,7),abs(camA_in_row(:,8)-height-1));
    water = water_A;

    % TONGUE
    ttongue_A = camA(:,10:11,:);
    lik_tongue_A = repmat(camA(:,12,:),[1, 2, 1]);
    tongue_A_flip_h = abs(ttongue_A(:,2,:)-height-1);
    tongue_A=cat(2,ttongue_A(:,1,:),tongue_A_flip_h);
    tongue_A(lik_tongue_A<likelihood_threshold) = nan;
    tongueA_in_row = cat(2,camA_in_row(:,10),abs(camA_in_row(:,11)-height-1));
    tongue = tongue_A;


    % Left or Right paw??
    % Right
    if front_tracked_flag==true
        tpawR_B = camB(:,1:2,:);
        lik_pawR_B = repmat(camB(:,3,:),[1, 2, 1]);
        pawR_B_flip_h = abs(tpawR_B(:,2,:)-heightB-1);
        pawR_B=cat(2,tpawR_B(:,1,:),pawR_B_flip_h);
        pawR_B(lik_pawR_B<likelihood_threshold) = nan;
        %Left
        tpawL_B = camB(:,4:5,:);
        lik_pawL_B = repmat(camB(:,6,:),[1, 2, 1]);
        pawL_B_flip_h = abs(tpawL_B(:,2,:)-heightB-1);
        pawL_B=cat(2,tpawL_B(:,1,:),pawL_B_flip_h);
        pawL_B(lik_pawL_B<likelihood_threshold) = nan;
        % Water
        twater_B = camB(:,7:8,:);
        lik_water_B = repmat(camB(:,9,:),[1, 2, 1]);
        wat_B_flip_h = abs(twater_B(:,2,:)-heightB-1);
        wat_B=cat(2,twater_B(:,1,:),pawL_B_flip_h);
        wat_B(lik_water_B<likelihood_threshold) = nan;
        %water_posB_px = [562, 413];
        %water_posB = [water_posB_px(1), abs(water_posB_px(2)-heightB-1)];

        water_posB = squeeze(round(median(wat_B,1,'omitnan')));
        if size(water_posB,1)==1, water_posB=water_posB'; end

        % Save
        save(strcat(out_mat_folder,filesep,'session_raw_data.mat'),'resnet','resnet_front', ...
            'width','height','widthB','heightB','likelihood_threshold','frame_rate', ...
            'csvs', 'n_trials', 'siz', 'max_siz', 'camA', 'camB',...
            'paw_A','pawA_in_row','snout_A','snoutA_in_row','water_A','waterA_in_row','water',...
            'tongue_A','tongueA_in_row',...
            'pawR_B','pawL_B','wat_B','water_posB','n_max_points',...
            '-v7.3');
    end

    figure
    if size(paw_A,3)==1, tt=1; else, tt=6; end
    subplot(121)
    plot(paw_A(1:siz(tt),1,tt),paw_A(1:siz(tt),2,tt),'-'); hold on
    plot(water_A(1:siz(tt),1,tt),water_A(1:siz(tt),2,tt),'.');
    plot(snout_A(1:siz(tt),1,tt),snout_A(1:siz(tt),2,tt),'.');
    plot(tongue_A(1:siz(tt),1,tt),tongue_A(1:siz(tt),2,tt),'.'); hold off
    legend('paw','water','snout','tongue','location','southeast','box','off')
    axis square
    if front_tracked_flag==true
        subplot(122)
        plot(pawR_B(1:siz(tt),1,tt),pawR_B(1:siz(tt),2,tt),'.'); hold on
        plot(pawL_B(1:siz(tt),1,tt),pawL_B(1:siz(tt),2,tt),'.');
        plot(wat_B(1:siz(tt),1,tt),wat_B(1:siz(tt),2,tt),'.','MarkerSize',10);
        plot(water_posB(1,tt),water_posB(2,tt),'*','MarkerSize',10); hold off
        legend('pawR','pawL','water','location','southeast','box','off')
        axis square
    end

else
    load(strcat(out_mat_folder,filesep,'session_raw_data.mat'));
    disp('loaded tracked features data')

end

%% Remove DLC jumps

paw = zeros(size(paw_A));
paw_sess = [];
thres_jmp = 65;

for t=1:n_trials
    % remove jumps
    paw(1:siz(t),1,t) = remove_DLCjumps(paw_A(1:siz(t),1,t),thres_jmp,0);
    %paw(1:siz(t),2,t) = remove_DLCjumps(paw_A(1:siz(t),2,t),thres_jmp,0);
    paw(1:siz(t),2,t) = paw_A(1:siz(t),2,t);

    % Join every trial in a row
    paw_sess = cat(1,paw_sess,...
        squeeze(paw(1:siz(t),:,t)));
end

figure
subplot(211)
plot(paw_sess(:,1));
subplot(212)
plot(paw_sess(:,2));


% Remove jumps and get speed
speed_paw = cat(1,nan(1,size(paw,2),size(paw,3)),...
    diff(paw));
N_v = 1;
K = 3;
F = 7;
speed_paw_s = nan(size(paw));
for t = 1:n_trials
    for dim = 1:2
        x = squeeze(paw(:,dim,t));
        speed_paw_s(:,dim,t) = applySGolayDerivationFilter(x, N_v, K, F);
    end
end

% 2D speed
dt = 1/frame_rate;
speed_partial_dt = speed_paw./dt;
speed_2D = squeeze(sqrt(...
    speed_partial_dt(:,1,:).^2 + ...
    speed_partial_dt(:,2,:).^2  ...
    ));

if size(paw_A,3)==1, tt=1; else, tt=6; end
figure
subplot(321)
plot(paw(1:siz(tt),1,tt));

subplot(322)
plot(speed_paw(1:siz(tt),1,tt)); hold on
plot(speed_paw_s(1:siz(tt),1,tt));

subplot(323)
plot(paw(1:siz(tt),2,tt));

subplot(324)
plot(speed_paw(1:siz(tt),2,tt)); hold on
plot(speed_paw_s(1:siz(tt),2,tt));

subplot(325)
plot(speed_2D(:,tt)); %xlim([0 180])

%% ----------------------------------------------------

%%  GO TRIAL BY TRIAL, FIND REACHES
warning('off','all')
close all

% Allocate space and variables
% Trial variables
nreaches_each_trial=zeros(n_trials,1); % n_reaches
trial_first_attempt=zeros(n_trials,1); %  first attempt
last_frame_wWater = zeros(n_trials,1);
time_to_hit = nan(n_trials,1);

% Params for find peaks (pk height defined by proximity to water)
% water position
water_x=water_A(:,1,:);
med_watx = nanmedian(water_x(:));
% Findpeaks params
prec_close_water = 0.75;
pkheight = floor(prec_close_water*med_watx); %288; % 299; %299
pkdist= 15; % 60;
pkprominence= 15 ; %25

% Find peaks across session to find IRI (interreach interval)
[pks_sess, pks_frames] = findpeaks(paw_sess(:,1),'MinPeakHeight',...
    pkheight,'MinPeakDistance',pkdist,'MinPeakProminence',pkprominence,...
    'Annotate','extents');

% parameters for defining if water was available (purposeful reaches)
extra_frames=10;
nan_min = 30;

% REACH PARAMETERS
% parameters of reach interval
reach_goBack_frames=40;
reach_goFront_frames=60;
max_reach=reach_goBack_frames+1;
time_range = ((0:(reach_goBack_frames+reach_goFront_frames))-reach_goBack_frames)/frame_rate;

% reach variables initialization
r_num=0;
reach_trial = [];
reach_range_mat = [];
reach_mat = [];
water_reach_mat = [];
reach_speed_mat = [];

r_paw_used_mat = [];

start_forw_mat = [];
no_turn_from_forw_mat = [];
dur_forw_mat = [];
dist_trav_forw_mat= [];
displacement_forw_mat = [];
tortuosity_forw_mat = [];

stop_reach_mat = [];
dur_start_stop_mat = [];
dist_trav_start_stop_mat = [];
displacement_start_stop_mat = [];
tortuosity_start_stop_mat = [];

median_speed_forward_mov_mat = [];
median_speed_forward_mov_start_mat = [];
max_speed_forward_mov_mat = [];

dist_trav_full_reach_mat = [];
displacement_full_reach_mat = [];
tortuosity_full_reach_mat = [];

cat_reach = [];
purpose_reach = [];
hit_reach = [];
success_reach = [];

% ------------------------------------------------------------------
% TRIAL LOOP
disp('Finding reaches...')

%tt=1
for tt = 1:n_trials
    disp(strcat('trial ',num2str(tt)))

    % WATER ------------- (consider improving at some point)
    % water availability
    t_water_x = water_A(1:siz(tt),1,tt);
    cumsum_nan = cumsum(isnan(t_water_x));

    if max(cumsum_nan) == numel(cumsum_nan) % only nans
        last_water = 1;
    elseif max(cumsum_nan)<nan_min
        last_water = find(cumsum_nan==max(cumsum_nan),1);
    elseif cumsum_nan(end) == nan_min
        last_water = siz(tt)-extra_frames;
    else
        last_water = find(cumsum_nan>nan_min,1);
    end
    purpose_length=1:(last_water+extra_frames);
    water_length=1:last_water;
    last_frame_wWater(tt,1)=last_water+extra_frames;
    time_to_hit(tt,1)=numel(water_length)/frame_rate;
    %figure, plot(cumsum_nan)
    % ---------------------------

    % CHECK PAW USED FOR REACHING
    if front_tracked_flag==true
        pawR_B_tt = pawR_B(1:siz(tt),:,tt);
        pawL_B_tt = pawL_B(1:siz(tt),:,tt);
        was_pos = repmat(squeeze(water_posB(:,tt))',[siz(tt) 1]);
        %was_pos=wat_B(:,:,tt);

        if show_plot==true
            figure
            subplot(211)
            plot(pawR_B_tt(:,1)); hold on
            plot(pawL_B_tt(:,1));
            plot(was_pos(:,1),'.'); hold off
            xlabel('frames'); ylabel('y')

            subplot(212)
            plot(pawR_B_tt(:,2)); hold on
            plot(pawL_B_tt(:,2));
            plot(was_pos(:,2),'.'); hold off
            xlabel('frames'); ylabel('z')
            legend('R','L','w','location','southeast')

        end
    end

    % DOMINANT PAW ----------------------------
    % trial paw position
    t_domPaw = paw(1:siz(tt),:,tt);
    t_domPaw_x = squeeze(t_domPaw(:,1));
    % trial paw speed
    t_speed_domPaw = speed_paw(1:siz(tt),:,tt);
    t_speed_domPaw_x = squeeze(t_speed_domPaw(:,1));
    t_speed_domPaw_x_s= squeeze(speed_paw_s(1:siz(tt),1,tt));
    % FIND PEAKS
    % findpeaks params adjustment to trials
    if pkdist>=siz(tt)
        pkdist=siz(tt)-2;
    end
    % find peaks
    [xmax, xmax_frame,w,p] = findpeaks(t_domPaw_x,'MinPeakHeight',...
        pkheight,'MinPeakDistance',pkdist,'MinPeakProminence',pkprominence);

    if show_plot==true
        figure
        subplot(211)
        findpeaks(t_domPaw_x,'MinPeakHeight',...
            pkheight,'MinPeakDistance',pkdist,'MinPeakProminence',pkprominence,...
            'Annotate','extents');
        hold on;
        plot(t_water_x,'k');
        plot(purpose_length,ones(1,numel(purpose_length)).*pkheight,'k');
        %plot(trial_nondomPaw_x(trial_length))
        %xline(xmax_frame); hold off
        xlim([0 siz(tt)]); ylim([20, 380]);
        title('position'); ylabel('x (px)'); xlabel('frames')
        subplot(212)
        plot(t_speed_domPaw_x); hold on
        plot(t_speed_domPaw_x_s,'color',[.5 .5 .5 0.6]); hold on
        plot(xmax_frame,t_speed_domPaw_x(xmax_frame),'*','markersize',2,'linewidth',5)
        %xline(xmax_frame); hold off
        xlim([0 siz(tt)]); grid on
        title('speed'); ylabel('x (px/frame)'); xlabel('frames')
        set(gcf, 'Position', [2182 251 1318 700])
        %if save_flag, saveas(gcf,strcat(folders.out,filesep,'trial_reaches','trial',num2str(tt),'.png'),'png'); end
    end


    % Trial params
    % n reaches
    nreaches_trial = numel(xmax);
    % first attempt?
    if nreaches_trial == 0
        isfirst_attempt = nan;
    elseif (nreaches_trial == 1 || xmax_frame(2)>purpose_length(end))
        isfirst_attempt = 1;
    else
        isfirst_attempt = 0;
    end
    trial_first_attempt(tt,1) = isfirst_attempt;
    nreaches_each_trial(tt,1) = nreaches_trial;




    %% LOOP OF REACHES
    close all
    if nreaches_trial == 0
        disp(strcat('no reaches in trial',num2str(tt)))
    else
        for r_ind = 1:nreaches_trial
            %r_ind=1
            % reach max frame
            rr=xmax_frame(r_ind);
            nn_max_points_tt = siz(tt);

            % grab frames from previous trial if reach was right at the beggiging of current trial
            if rr<=reach_goBack_frames
                last_frame_prev_tt1 = siz(tt-1);
                if last_frame_prev_tt1<reach_goBack_frames-rr
                    last_frame_prev_tt2 = siz(tt-2);
                    last_trials=cat(1,paw(1:last_frame_prev_tt2,:,tt-2),paw(1:last_frame_prev_tt1,:,tt-1));
                    last_trials_speed=cat(1,speed_paw(1:last_frame_prev_tt2,:,tt-2),speed_paw(1:last_frame_prev_tt1,:,tt-1));
                    last_trials_water=cat(1,water_A(1:last_frame_prev_tt2,:,tt-2),water_A(1:last_frame_prev_tt1,:,tt-1));
                    last_trials_pawR_B=cat(1,pawR_B(1:last_frame_prev_tt2,:,tt-2),pawR_B(1:last_frame_prev_tt1,:,tt-1));
                    last_trials_pawL_B=cat(1,pawL_B(1:last_frame_prev_tt2,:,tt-2),pawL_B(1:last_frame_prev_tt1,:,tt-1));
                else
                    last_trials=paw(1:last_frame_prev_tt1,:,tt-1);
                    last_trials_speed=speed_paw(1:last_frame_prev_tt1,:,tt-1);
                    last_trials_water=water_A(1:last_frame_prev_tt1,:,tt-1);
                    if front_tracked_flag==true
                        last_trials_pawR_B=pawR_B(1:last_frame_prev_tt1,:,tt-1);
                        last_trials_pawL_B=pawL_B(1:last_frame_prev_tt1,:,tt-1);
                    end
                end
                last_siz=size(last_trials_water,1);
                reach_prev = last_siz-(reach_goBack_frames-rr):last_siz;
                reach_curr = 1:rr+reach_goFront_frames;
                reach_range=nan(1,reach_goBack_frames+reach_goFront_frames+1);
                reach_range(reach_goBack_frames-rr+2:end)=reach_curr;
                reach=cat(1,squeeze(last_trials(reach_prev,:)),squeeze(paw(reach_curr,:,tt)));
                reach_speed=cat(1,squeeze(last_trials_speed(reach_prev,:)),squeeze(speed_paw(reach_curr,:,tt)));
                water_reach=cat(1,squeeze(last_trials_water(reach_prev,:)),squeeze(water_A(reach_curr,:,tt)));
                if front_tracked_flag==true
                    pawR_B_reach = cat(1,squeeze(last_trials_pawR_B(reach_prev,:)),squeeze(pawR_B(reach_curr,:,tt)));
                    pawL_B_reach = cat(1,squeeze(last_trials_pawL_B(reach_prev,:)),squeeze(pawL_B(reach_curr,:,tt)));
                end

            elseif nn_max_points_tt-rr<reach_goFront_frames
                reach_curr = rr-reach_goBack_frames:nn_max_points_tt;
                reach_nex = 1:reach_goFront_frames-(nn_max_points_tt-rr);
                reach_range=nan(1,reach_goBack_frames+reach_goFront_frames+1);
                reach_range(1:numel(reach_curr))=reach_curr;
                reach_range(numel(reach_curr)+1:end)=reach_nex;
                reach=cat(1,squeeze(paw(reach_curr,:,tt)),squeeze(paw(reach_nex,:,tt+1)));
                reach_speed=cat(1,squeeze(speed_paw(reach_curr,:,tt)),squeeze(speed_paw(reach_nex,:,tt+1)));
                water_reach=cat(1,squeeze(paw(reach_curr,:,tt)),squeeze(paw(reach_nex,:,tt+1)));
                if front_tracked_flag==true
                    pawR_B_reach=cat(1,squeeze(pawR_B(reach_curr,:,tt)),squeeze(pawR_B(reach_nex,:,tt+1)));
                    pawL_B_reach=cat(1,squeeze(pawL_B(reach_curr,:,tt)),squeeze(pawL_B(reach_nex,:,tt+1)));
                end
            else
                reach_range=rr-reach_goBack_frames:rr+reach_goFront_frames;
                reach_curr = reach_range;
                reach=squeeze(paw(reach_range,:,tt));
                reach_speed=squeeze(speed_paw(reach_range,:,tt));
                water_reach=squeeze(water_A(reach_range,:,tt));
                if front_tracked_flag==true
                    pawR_B_reach=squeeze(pawR_B(reach_range,:,tt));
                    pawL_B_reach=squeeze(pawL_B(reach_range,:,tt));
                end
            end


            %% SPEED and DISTANCE TRAVELLED
            % speed
            speed_partial_dt = reach_speed./dt;
            speed_2D_reach = squeeze(sqrt(...
                speed_partial_dt(:,1,:).^2 + ...
                speed_partial_dt(:,2,:).^2  ...
                ));

            % distance travelled
            dist_points = squeeze(sqrt(...
                reach_speed(:,1,:).^2 + ...
                reach_speed(:,2,:).^2 ...
                ));
            dist_trav_reach = cumsum(dist_points,'omitnan');


            if show_plot==true
                clrs = lines(5);
                figure
                subplot(411)
                plot(t_domPaw_x,'k'); hold on
                plot(t_water_x,'color',[.5 .5 .5]);
                hold on
                plot(reach_curr,reach(1:numel(reach_curr),1),'color',clrs(1,:),'linewidth',1.5);
                plot(xmax_frame,xmax,'*','LineWidth',5,'color',clrs(2,:))
                plot(purpose_length,ones(1,numel(purpose_length)).*pkheight,'k');
                hold off
                ylabel('x')
                legend('dom paw','peak','proeminence','width','water_A')

                subplot(412)
                plot(t_speed_domPaw_x,'k','linewidth',1)
                hold on
                plot(reach_curr,reach_speed(1:numel(reach_curr),1),'color',clrs(1,:),'linewidth',1.5)
                plot(xmax_frame,t_speed_domPaw_x(xmax_frame),'*','LineWidth',5,'color',clrs(2,:))
                hold off
                legend('trial','peaks','selected reach');
                ylabel('speed x')

                subplot(413)
                plot(t_domPaw(:,2),'k','linewidth',1)
                hold on
                plot(reach_curr,reach(1:numel(reach_curr),2),'color',clrs(1,:),'linewidth',1.5)
                plot(xmax_frame,t_domPaw(xmax_frame,2),'*','LineWidth',5,'color',clrs(2,:))
                hold off
                ylabel('z')
                set(gcf, 'Position', [2262 380 1367 536])


                subplot(413)
                plot(t_domPaw(:,2),'k','linewidth',1)
                hold on
                plot(reach_curr,reach(1:numel(reach_curr),2),'color',clrs(1,:),'linewidth',1.5)
                plot(xmax_frame,t_domPaw(xmax_frame,2),'*','LineWidth',5,'color',clrs(2,:))
                hold off
                ylabel('z')

                if front_tracked_flag==true
                    subplot(414)
                    plot(pawR_B_tt(:,1)); hold on
                    plot(pawL_B_tt(:,1));
                    plot(was_pos(:,1),'o');
                    plot(reach_curr,pawR_B_reach(1:numel(reach_curr),1),'color',clrs(1,:),'linewidth',2)
                    plot(reach_curr,pawL_B_reach(1:numel(reach_curr),1),'color',clrs(2,:),'linewidth',2)
                    plot(xmax_frame,was_pos(xmax_frame,1),'*','LineWidth',5,'color',clrs(3,:))
                    legend('R','L','w')
                    ylabel('front cam')
                end
                set(gcf, 'Position', [2262  110 1367 861])
                %if save_fig_flag, saveas(gcf,strcat(folders.out,filesep,'detect_reach',num2str(r_ind),'trial',num2str(tt),'.png'),'png'); end
            end

            % CHECK WHICH PAW WAS USED --------------------------
            if front_tracked_flag==true
                check_paw_frames_before = 8;
                check_paw_frames_after = 5;
                check_paw_range = max_reach-check_paw_frames_before:max_reach+check_paw_frames_after;
                max_dist = 100;

                med_position_R = median(pawR_B_reach(check_paw_range,1),'omitnan');
                med_position_L = median(pawL_B_reach(check_paw_range,1),'omitnan');

                distance_RW = abs(med_position_R-(was_pos(1,1)));
                distance_LW = abs(med_position_L-(was_pos(1,1)));

                if (distance_RW<distance_LW || (isnan(distance_LW) && distance_RW<max_dist))
                    reach_paw_used = 'R';
                    r_paw_used = 1;
                elseif (distance_RW>distance_LW || (isnan(distance_RW) && distance_LW<max_dist))
                    reach_paw_used = 'L';
                    r_paw_used = 0;
                else
                    reach_paw_used = 'unknown';
                    r_paw_used = nan;
                end

                if show_plot==true
                    figure
                    plot(pawR_B_reach(check_paw_range,1),'color',clrs(1,:),'linewidth',2); hold on
                    plot(repmat(med_position_R,[length(check_paw_range) 1]),'color',clrs(1,:));
                    plot(pawL_B_reach(check_paw_range,1),'color',clrs(2,:),'linewidth',2)
                    plot(repmat(med_position_L,[length(check_paw_range) 1]),'color',clrs(2,:));
                    plot(was_pos(check_paw_range,1),'color',clrs(3,:),'linewidth',2)
                    plot(check_paw_frames_before+1,was_pos(check_paw_frames_before+1),'*','LineWidth',5,'color',clrs(3,:))
                    title(sprintf('%s%s','Used paw: ',reach_paw_used));
                    set(gcf,'position', [3183 375 560 420])
                end
            else
                reach_paw_used='?';
            end

            %-----------------------------------------------------
            % ONSET OF REACH: FORWARD MOVEMENT & FROM STOP
            % FORWARD MOVEMENT ISOLATED
            % Check if paw started in resting bar/inverted the direction of reach
            min_reach_time = 5; % 60-70 ms
            reaching_frames_minimum = reach_goBack_frames-min_reach_time;
            min_frames_down = 3; % 10-15 ms
            close_zero = 1;
            % Find foward movement start
            t2 = reach_speed(:,1);
            zt2=zeros(size(t2)); zt2(t2>close_zero)=1; zt2(reaching_frames_minimum:end)=1;
            fzt2 = flip(zt2,1); fzt2(end)=1; % Flip zt2 to find the first stable rise timepoint
            if isempty(find(diff(fzt2)~=0,1)), fzt2(end)=0; fzt2(end-2)=0; end
            fjumps_loc=find(diff(fzt2)~=0); fjumps_down_loc=fjumps_loc((1:2:end)); % find the position of jumps down
            diff_jumps=diff(fjumps_loc); inv_rise_jumps_diff = diff_jumps(1:2:end); % find number of frames after jump
            inv_stable_rising_start=fjumps_down_loc(inv_rise_jumps_diff>min_frames_down); % Impose a min number of frames to consider stalbe rise
            if isempty(inv_stable_rising_start)
                [max_inv_steps,max_inv_steps_loc]=max(inv_rise_jumps_diff);
                inv_stable_rising_start=fjumps_down_loc(max_inv_steps_loc);
            end
            rising_start=numel(fzt2)-inv_stable_rising_start(1); % Invert back
            rising_start = rising_start + 1; % correct for the diff shift

            % Paw turns back?
            if numel(fzt2)-fjumps_down_loc(1)==rising_start
                no_turns_from_forw = 1;
            else
                no_turns_from_forw = 0;
            end
            % Duration of the continuously forward movement
            dur_forw = numel(rising_start:max_reach)/frame_rate;
            % Path length (distance travelled), displacement and tortuosity
            [dist_trav_forw, displacement_forw, tortuosity_forw] = tortuosity_calc_2D(reach,[rising_start,max_reach]);

            % Max and media speed
            median_speed_forward_mov = median(reach_speed(rising_start:max_reach,1),'omitnan');
            median_speed_forward_mov_start = median(reach_speed(rising_start+5:max_reach-10,1),'omitnan');
            max_speed_forward_mov = max(reach_speed(rising_start:max_reach,1));



            % --------------------------------------------------------
            % REACH ENDING - STOP  -----------------------------------------
            close_zero_stop = 2;
            zt2_stop=zeros(size(t2)); zt2_stop(t2>close_zero_stop)=1; zt2_stop(1:max_reach)=0;
            stop_reach = find(zt2_stop==1,1);
            if isempty(stop_reach)
                stop_reach = max_reach+find(isnan(t2(max_reach:end)),1);
                stop_reach = stop_reach-1;
                if isempty(stop_reach)
                    stop_reach = length(t2);
                end
            end
            % Duration and displacement from rising to stop reach
            dur_start_stop = numel(rising_start:stop_reach)/frame_rate;
            duration_up_paw =  numel(max_reach:stop_reach)/frame_rate;
            [dist_trav_start_stop, displacement_start_stop, tortuosity_start_stop] = tortuosity_calc_2D(reach,[rising_start,stop_reach]);
            median_position_reach = median(reach(rising_start:stop_reach,1),'omitnan');
            median_post_reach = median(reach(max_reach:stop_reach,1),'omitnan');
            max_position_y_reach = reach(max_reach,2);

            % Distance, displacement and tortuosity of full reach
            [dist_trav_full_reach, displacement_full_reach, tortuosity_full_reach] = tortuosity_calc_2D(reach,...
                [find(~isnan(reach(:,1)),1),find(~isnan(reach(:,1)),1,'last')]);

            if show_plot==true
                figure
                subplot(411)
                plot(reach(:,1));
                xline(max_reach);
                xline(rising_start);
                xline(stop_reach,'g');
                xlim([0 100]);
                title('position in X')
                subplot(412)
                plot(reach_speed(:,1)); hold on
                plot(t2);
                plot(zt2,'.')
                xline(rising_start);
                xline(stop_reach,'g');
                xline(max_reach);  hold off
                xlim([0 100]);
                title('speed in X, forward mov start')
                subplot(413)
                plot(reach_speed(:,1)); hold on
                xline(max_reach);
                xline(stop_reach,'g');
                xlim([0 100]); hold off
                title('speed in X, movement in x start')
                subplot(414)
                plot(speed_2D_reach); hold on
                xline(max_reach);
                xline(stop_reach,'g');
                xlim([0 100]);
                title('speed 3D, overall movement start')
                set(gcf,'position',[680   160   777   818])
            end


            %% CATEGORIES
            % Reach Category
            % 1 -> balistic
            % 2 -> holding spout
            % 3 -> hiden start

            % Hidden start
            rising_start_min_pos = 125;
            % Duration for not balistic
            %dur_max_balist = 0.3;
            % Median position min to be considered holding
            med_position_min_for_hold = 250;

            if reach(rising_start,1) > rising_start_min_pos % Has hidden start
                which_cat = 3; % hidden start
            elseif median_post_reach > med_position_min_for_hold
                which_cat = 2; % holding spout
            else
                which_cat = 1; % balistic
            end


            %% Reach Class
            % Threshold outside ROI - successful?
            threshold_x_max = 300; threshold_x_min = 100;
            threshold_z_max = 270; threshold_z_min = 200; %150
            max_outliers=0;

            % PURPOSFULL OR NOT?
            if rr < purpose_length(end)
                is_purpose_reach = 1;

                % HIT OR MISS?
                if ( (rr == xmax_frame(nreaches_trial) && purpose_length(end)<n_max_points) ...
                        || (rr ~= xmax_frame(nreaches_trial) && xmax_frame(r_ind+1)>purpose_length(end)) )
                    is_hit = 1;

                    % SUCCESS OR FAILURE?
                    % Range around reach
                    if rr<=50
                        aroud_reach=1:rr+50;
                    elseif rr+100>n_max_points
                        aroud_reach=rr-50:n_max_points;
                    else
                        aroud_reach=rr-50:rr+100;
                    end

                    % Plot
                    if show_plot==true
                        figure(8)
                        subplot(211)
                        plot(water(:,1,tt),'.'); ylabel('x')
                        hold on
                        %plot(xmax_frame(r_ind),median_water(1),'o','LineWidth',2)
                        %plot(aroud_reach(1)-1+pos_out_x,water(aroud_reach(1)-1+pos_out_x,1,tt),'.')
                        yline(threshold_x_max,'g'); yline(threshold_x_min,'g');
                        hold off

                        subplot(212)
                        plot(water(:,2,tt),'.'); ylabel('z')
                        yline(threshold_z_max,'g'); yline(threshold_z_min,'g');

                        %saveas(gcf,strcat(folders.out,filesep,'successORmiss_trial',num2str(tt),'.png'),'png');
                    end

                    w_ar=squeeze(water(aroud_reach,:,tt));
                    z_w_ar=w_ar(w_ar(:,2)<threshold_z_min,2);
                    x_w_ar=w_ar(w_ar(:,1)>threshold_x_max,1);
                    if (...
                            numel(find(w_ar(:,1)>threshold_x_max))>0 || numel(find(w_ar(:,1)<threshold_x_min))>0 || ...  % frames with water outside ROI or to the sides
                            numel(find(w_ar(:,2)>threshold_z_max))>0 || (numel(find(w_ar(:,2)<threshold_z_min))>0 && min(z_w_ar)==z_w_ar(end)) )
                        %    (max(count_outROI)>=max_outliers && numel(find(creount_outROI==0))<=1) || ...
                        is_success = 0;
                    else
                        is_success = 1;
                    end
                else
                    is_hit = 0;
                    is_success = nan;
                end
            else
                is_purpose_reach = 0;
                is_hit = nan;
                is_success = nan;
            end


            %% FIGURE: REACH
            % plot
            if show_overview_reach_plot==true
                close(figure(10))
                fig_reach_prop = figure(10);
                set(fig_reach_prop,figOpt{:})
                speed_lim = [-2000 3000];
                speed_lim2 = [-100 5000];



                annotation('textbox', [0.07, 0.78, 0.15, 0.16],'String',{...
                    "TRIAL " + tt + " , REACH " + r_ind,...
                    "Paw used: " + reach_paw_used,...
                    }, 'FitBoxToText','on', 'EdgeColor', 'w');


                if which_cat==3
                    annotation('textbox', [0.07, 0.76, 0.75, 0.11],'String',[...
                        "CATEGORY:", "- Hidden start (3)",...
                        ], 'FitBoxToText','on', 'EdgeColor', 'w');
                elseif which_cat==1
                    annotation('textbox', [0.07, 0.76, 0.75, 0.11],'String',[...
                        "CATEGORY:", "- balistic reaching (1)",...
                        ], 'FitBoxToText','on', 'EdgeColor', 'w');
                elseif which_cat==2
                    annotation('textbox',[0.07, 0.76, 0.75, 0.11],'String',[...
                        "CATEGORY:", "- holding spout reaching (2)",...
                        ], 'FitBoxToText','on', 'EdgeColor', 'w');
                end


                annotation('textbox', [0.07, 0.63, 0.15, 0.16],'String',{...
                    "POSITION AND DURATIONS:",...
                    "Med x paw position after reach: " + sprintf('%.0f',median_position_reach) + " px (thres = "+ med_position_min_for_hold + ")"...
                    "Max z paw position at reach: " + sprintf('%.0f',max_position_y_reach)+ " px",...
                    "",...
                    "Reach foward duration: " + sprintf('%.3f',dur_forw) + " s",...
                    "Reach full duration: " + sprintf('%.3f',dur_start_stop) + " s",...
                    "",...
                    "Med speed forward reach: " + sprintf('%.3f',median_speed_forward_mov) + " px/s", ...
                    "Med speed forward reach start: " + sprintf('%.3f',median_speed_forward_mov_start) + " px/s", ...
                    "Max speed forward reach: " + sprintf('%.3f',max_speed_forward_mov) + " px/s", ...
                    }, 'FitBoxToText','on', 'EdgeColor', 'w');


                annotation('textbox', [0.07, 0.37, 0.15, 0.16],'String',{...
                    "CLASSIFICATION",...
                    "purposeful (or not)? " + sprintf('%.0f',is_purpose_reach),...
                    "hit (or miss)? " + sprintf('%.0f',is_hit),...
                    "success (or fail)? " + sprintf('%.0f',is_success),"",...
                    }, 'FitBoxToText','on', 'EdgeColor', 'w');


                annotation('textbox', [0.07, 0.28 ...
                    , 0.75, 0.11],'String',{...
                    "FORWARD BALISTIC MOVEMENT",...
                    "rising forward = " + sprintf('%.3f',time_range(rising_start)),...
                    "no turns from forw = " + no_turns_from_forw,...
                    "distance trav = " + sprintf('%.1f',dist_trav_forw) + " px",...
                    "displacement = " + sprintf('%.1f',displacement_forw)+ " px",...
                    "tortuosity = " + sprintf('%.3f',tortuosity_forw),...
                    "",...
                    "MOVEMENT FULL REACH",...
                    "stop time = " + sprintf('%.3f',time_range(stop_reach)),...
                    "distance trav start-stop = " +  sprintf('%.1f',dist_trav_start_stop) + " px",...
                    "distance trav = " + sprintf('%.1f',dist_trav_full_reach) + " px",...
                    "displacement = " + sprintf('%.1f',displacement_full_reach)+ " px",...
                    "tortuosity = " + sprintf('%.3f',tortuosity_full_reach),...
                    }, 'FitBoxToText','on', 'EdgeColor', 'w')



                subplot(4,3,2)
                plot(time_range,reach(:,1),'linewidth',1.5); hold on
                xline(time_range(max_reach),'--','color',clrs(3,:));
                xline(time_range(rising_start),'--','color',[.7 .7 .7]);
                xline(time_range(stop_reach),'--','color',[.7 .7 .7]);
                yline(rising_start_min_pos,'--','color',[.96 .95 .95]);
                ylabel('x (px)'); xlabel('time (sec)'); hold off
                title('Position')
                axis([min(time_range) max(time_range) 20 380])
                set(gca,axeOpt{:})

                subplot(4,3,3)
                plot(time_range,reach_speed(:,1)./dt,'linewidth',1.5); hold on
                xline(time_range(max_reach),'--','color',clrs(3,:));
                xline(time_range(rising_start),'--','color',[.7 .7 .7]);
                xline(time_range(stop_reach),'--','color',[.7 .7 .7]);
                yline(0,'--','color',[.9 .9 .9]);
                ylim(speed_lim)
                ylabel('x (px/s)'); xlabel('time (sec)'); hold off
                title('Speed')
                xlim([min(time_range) max(time_range)])
                set(gca,axeOpt{:})

                subplot(4,3,5)
                plot(time_range,reach(:,2),'linewidth',1.5); hold on
                xline(time_range(max_reach),'--','color',clrs(3,:));
                xline(time_range(rising_start),'--','color',[.7 .7 .7]);
                xline(time_range(stop_reach),'--','color',[.7 .7 .7]);
                ylabel('z (px)'); xlabel('time (sec)'); hold off
                axis([min(time_range) max(time_range) 70 400])
                set(gca,axeOpt{:})

                subplot(4,3,6)
                plot(time_range,reach_speed(:,2)./dt,'linewidth',1.5); hold on
                yline(0,'--','color',[.9 .9 .9]);
                xline(time_range(max_reach),'--','color',clrs(3,:));
                xline(time_range(rising_start),'--','color',[.7 .7 .7]);
                xline(time_range(stop_reach),'--','color',[.7 .7 .7]);
                ylabel('z (px/s)'); xlabel('time (sec)'); hold off
                xlim([min(time_range) max(time_range)])
                ylim(speed_lim)
                set(gca,axeOpt{:})

                subplot(4,3,8)
                plot(time_range,dist_trav_reach,'Color',clrs(1,:),'linewidth',1.5)
                xline(time_range(max_reach),'--','color',clrs(3,:));
                xline(time_range(rising_start),'--','color',[.7 .7 .7]);
                xline(time_range(stop_reach),'--','color',[.7 .7 .7]);
                ylabel('distance travelled (px)'); xlabel('time (sec)'); hold off
                xlim([min(time_range) max(time_range)])
                set(gca,axeOpt{:})
                title('Distance Travelled')

                subplot(4,3,9)
                plot(time_range,speed_2D_reach,'Color',clrs(1,:),'linewidth',1.5)
                xline(time_range(max_reach),'--','color',clrs(3,:));
                xline(time_range(rising_start),'--','color',[.7 .7 .7]);
                xline(time_range(stop_reach),'--','color',[.7 .7 .7]);
                ylabel('speed (px/s)'); xlabel('time (sec)'); hold off
                xlim([min(time_range) max(time_range)])
                ylim(speed_lim2)
                set(gca,axeOpt{:})
                title('Absolut Speed 2D')

                subplot(4,3,[11,12])
                plot(t_domPaw_x,'color', [.5 .5 .5]); hold on
                plot(t_water_x,'color','k');
                hold on
                plot(reach_curr,reach(1:numel(reach_curr),1),'color',clrs(1,:),'linewidth',1.5);
                plot(xmax_frame,xmax,'*','LineWidth',1,'color',clrs(3,:))
                plot(purpose_length,ones(1,numel(purpose_length)).*pkheight,'--k');
                set(gca,axeOpt{:})
                xlabel('frames in trial'); ylabel('x (px)')
                title('Reach position in trial')

                set(fig_reach_prop,'Position',[2124 89 1124 779])
                if save_fig_flag, saveas(gcf, strcat(output_folder,filesep,'reach',num2str(r_ind),'_trial',num2str(tt),'.png'),'png'); end

            end


            %% Update Reach Variables
            r_num=r_num+1;
            reach_trial = cat(1,reach_trial,tt);
            reach_range_mat = cat(1,reach_range_mat,reach_range);
            reach_mat = cat(3,reach_mat,reach);
            water_reach_mat = cat(3,water_reach_mat,water_reach);
            reach_speed_mat = cat(3,reach_speed_mat,reach_speed);
            if front_tracked_flag==true
                r_paw_used_mat = cat(1,r_paw_used_mat,r_paw_used);
            end

            start_forw_mat = cat(1,start_forw_mat,rising_start);
            no_turn_from_forw_mat = cat(1,no_turn_from_forw_mat,no_turns_from_forw);
            dur_forw_mat = cat(1,dur_forw_mat,dur_forw);
            dist_trav_forw_mat = cat(1,dist_trav_forw_mat,dist_trav_forw);
            displacement_forw_mat = cat(1,displacement_forw_mat,displacement_forw);
            tortuosity_forw_mat = cat(1,tortuosity_forw_mat,tortuosity_forw);

            stop_reach_mat = cat(1,stop_reach_mat,stop_reach);
            dur_start_stop_mat = cat(1,dur_start_stop_mat,dur_start_stop);
            dist_trav_start_stop_mat = cat(1,dist_trav_start_stop_mat,dist_trav_start_stop);
            displacement_start_stop_mat = cat(1,displacement_start_stop_mat,displacement_start_stop);
            tortuosity_start_stop_mat = cat(1,tortuosity_start_stop_mat,tortuosity_start_stop);

            dist_trav_full_reach_mat = cat(1,dist_trav_full_reach_mat,dist_trav_full_reach);
            displacement_full_reach_mat = cat(1,displacement_full_reach_mat,displacement_full_reach);
            tortuosity_full_reach_mat = cat(1,tortuosity_full_reach_mat,tortuosity_full_reach);

            median_speed_forward_mov_mat = cat(1,median_speed_forward_mov_mat,median_speed_forward_mov);
            median_speed_forward_mov_start_mat = cat(1,median_speed_forward_mov_start_mat,median_speed_forward_mov_start);
            max_speed_forward_mov_mat = cat(1,max_speed_forward_mov_mat,max_speed_forward_mov);

            cat_reach = cat(1,cat_reach,which_cat);
            purpose_reach = cat(1,purpose_reach,is_purpose_reach);
            hit_reach = cat(1,hit_reach,is_hit);
            success_reach = cat(1,success_reach,is_success);



        end % end of if there are reaches in trial
    end % end of reach loop
end % end trial loop

%% SAVE VARIABLES

sessions.video.width = width;
sessions.video.height = height;
sessions.videoB.widthB = widthB;
sessions.videoB.heightB = heightB;
sessions.video.resnet = resnet;
sessions.videoB.resnetB = resnet_front;
sessions.video.frame_rate = frame_rate;
sessions.video.likelihood_threshold = likelihood_threshold;
% session data params
sessions.mouse_info.phenotype = phenotype;
sessions.params.interpol_clean_params.thres_jmp = thres_jmp;
sessions.params.derivative.N = N_v;
sessions.params.derivative.K = K;
sessions.params.derivative.F = F;

% Session data variables
sessions.session = session;
sessions.paw_ss = paw_sess;
sessions.pks_sess = pks_sess;
sessions.pks_frames = pks_frames;
sessions.dt = dt;

% trial data variables
trials.paw = paw;
trials.speed_paw = speed_paw;
trials.speed_paw_s = speed_paw_s;
trials.speed_partial_dt = speed_partial_dt;
trials.speed_2D = speed_2D;
trials.water = water;
trials.snout = snout_A;

% other trial varaibles
trials.nreaches_each_trial = nreaches_each_trial;
trials.first_attempt = trial_first_attempt;
trials.last_frame_wWater = last_frame_wWater;
trials.time_to_hit = time_to_hit;
trials.n_trials = n_trials;
trials.size_trial = siz;
% Trial paramers
trials.trial_params.find_peaks.prec_close_water = prec_close_water;
trials.trial_params.find_peaks.pkheight = pkheight;
trials.trial_params.find_peaks.pkdist = pkdist;
trials.trial_params.find_peaks.pkprominence = pkprominence;
trials.trial_params.water_availability.extra_frames = extra_frames;
trials.trial_params.water_availability.nan_min = nan_min;

% reaches data variables
reaches.r_num = r_num;
reaches.reach_trial =reach_trial;
reaches.reach_range_mat = reach_range_mat;
reaches.reach_mat = reach_mat;
reaches.water_reach_mat = water_reach_mat;
reaches.reach_speed_mat = reach_speed_mat;
reaches.time_range = time_range;

% start forw
reaches.start_forw.start_forw_mat = start_forw_mat;
reaches.start_forw.no_turn_from_forw_mat = no_turn_from_forw_mat;
reaches.start_forw.dur_forw_mat = dur_forw_mat;
reaches.start_forw.dist_trav_forw_mat = dist_trav_forw_mat;
reaches.start_forw.displacement_forw_mat = displacement_forw_mat;
reaches.start_forw.tortuosity_forw_mat = tortuosity_forw_mat;

% start to stop
reaches.start_stop.stop_reach_mat = stop_reach_mat;
reaches.start_stop.dur_start_stop_mat = dur_start_stop_mat;
reaches.start_stop.dist_trav_start_stop_mat = dist_trav_start_stop_mat;
reaches.start_stop.displacement_start_stop_mat = displacement_start_stop_mat;
reaches.start_stop.tortuosity_start_stop_mat = tortuosity_start_stop_mat;

% full reach
reaches.full_reach.dist_trav_full_reach_mat = dist_trav_full_reach_mat;
reaches.full_reach.displacement_full_reach_mat = displacement_full_reach_mat;
reaches.full_reach.tortuosity_full_reach_mat = tortuosity_full_reach_mat;

% med & max speed
reaches.speed_meds.median_speed_forward_mov_mat = median_speed_forward_mov_mat;
reaches.speed_meds.median_speed_forward_mov_start_mat = median_speed_forward_mov_start_mat;
reaches.speed_meds.max_speed_forward_mov_mat = max_speed_forward_mov_mat;

% reach category/class
reaches.cat_reach = cat_reach;
reaches.purpose_reach = purpose_reach;
reaches.hit_reach = hit_reach;
reaches.success_reach = success_reach;
reaches.Rpaw_used = r_paw_used_mat;

% reaches params
reaches.reach_params.reach_interval.reach_goBack_frames = reach_goBack_frames;
reaches.reach_params.reach_interval.reach_goFront_frames = reach_goFront_frames;
reaches.reach_params.reach_interval.max_reach = max_reach;
reaches.reach_params.start_reach.min_reach_time = min_reach_time;
reaches.reach_params.start_reach.reaching_frames_minimum = reaching_frames_minimum;
reaches.reach_params.start_reach.min_frames_down = min_frames_down;
reaches.reach_params.start_reach.close_zero = close_zero;
reaches.reach_params.start_reach.close_zero_stop  = close_zero_stop;

reaches.reach_params.cat_and_class.rising_start_min_pos = rising_start_min_pos;
reaches.reach_params.cat_and_class.med_position_min_for_hold = med_position_min_for_hold;
reaches.reach_params.cat_and_class.succ.threshold_x_max = threshold_x_max;
reaches.reach_params.cat_and_class.succ.threshold_x_min = threshold_x_min;
reaches.reach_params.cat_and_class.succ.threshold_z_max = threshold_z_max;
reaches.reach_params.cat_and_class.succ.threshold_z_min = threshold_z_min;
reaches.reach_params.cat_and_class.succ.max_outliers = max_outliers;




save(strcat(out_mat_folder,filesep,'session_reaching_data.mat'),'reaches','trials','sessions','mouse_info','-v7.3')
disp('Finished session!')


%end % close session loop
disp('All done!!!')









