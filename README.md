# post-collision-velocity-and-Velocity-ratio-analysis
for thesis
clc
clear;close all 
factor=5.1372e-6;       % 照片上一个像素代表的实际距离 7.2X
% factor=5.1372e-6*7.2/5; % 照片上一个像素代表的实际距离 4X
T=1/25000;              % 一个时间间隔
t=T;
Rho=750;                % 液滴的密度
sigma=0.025;            % 液滴的表面张力系数
interval=4;             % 人为选定照片张数
interval_post=4;        % number of frames after separation used for post-collision velocity
gap=0;                  % 为读取中间图片的情况预备
cropHeight = 672;       % 确定裁剪的高度
Highbegin = 2;
cropWidth = 768; 
Widthbegin = 1;
path='D:\rawam\research work\thesis\velocity\5 atm\bouncing 1\';

for k=10:1:10
    radii_left_all=[];
    radii_right_all=[];
    centers_all=[];
    We_all=[];
    B_all=[];
    relative_velocity_all=[];
    velocity_left=[];
    velocity_right=[];
    real_diameter_all=[];

    % ── post-collision storage ──────────────────────────────────────────
    centers_post=[];          % stores droplet centres after separation
    radii_post_left_all=[];
    radii_post_right_all=[];
    velocity_post_left_all=[];
    velocity_post_right_all=[];
    relative_velocity_post_all=[];
    % ────────────────────────────────────────────────────────────────────

    dirname=['time' num2str(k)];
    files = dir(fullfile(path,dirname,'\add1\','*.tif')); 
    lengthFiles = length(files);

    % ════════════════════════════════════════════════════════════════════
    % STEP 1 — find the frame where collision starts (two → one region)
    % ════════════════════════════════════════════════════════════════════
    collision_start = [];
    for i = 1:lengthFiles
        Img_before = imread(strcat(path,dirname,'\add1\',files(i+gap).name));
        Img = Img_before(Highbegin:cropHeight, Widthbegin:cropWidth);
        Img1 = bwareaopen(Img,500,8);
        Img2 = imcomplement(Img1);
        Img3 = bwareaopen(Img2,100,8);
        figure(k),imshow(Img3);
        [B,L] = bwboundaries(Img3,'holes');
        if length(B)==1       % droplets have merged → collision started
            collision_start = i;
            break
        end
    end
    number(k,:) = i-1;       % last pre-collision frame index

    if number(k,:)>4
        interval=4;
    else
        interval=number(k,:)-1;
    end

    % ════════════════════════════════════════════════════════════════════
    % STEP 2 — find the frame where separation occurs (one → two regions)
    %          i.e. liquid bridge breaks, droplets are free again
    % ════════════════════════════════════════════════════════════════════
    separation_frame = [];
    if ~isempty(collision_start)
        for i = collision_start:lengthFiles
            Img_before = imread(strcat(path,dirname,'\add1\',files(i+gap).name));
            Img = Img_before(Highbegin:cropHeight, Widthbegin:cropWidth);
            Img1 = bwareaopen(Img,500,8);
            Img2 = imcomplement(Img1);
            Img3 = bwareaopen(Img2,100,8);
            [B_sep,~] = bwboundaries(Img3,'holes');
            if length(B_sep) >= 2   % two separate regions → bridge broken
                separation_frame = i;
                fprintf('Separation detected at frame %d\n', separation_frame);
                break
            end
        end
    end

    % ════════════════════════════════════════════════════════════════════
    % STEP 3 — pre-collision: track droplets in frames 1..number(k)
    % ════════════════════════════════════════════════════════════════════
    for i=1:number(k)
        Img_origin_before_crop = imread(strcat(path,dirname,'\add1\',files(i+gap).name));
        Img_origin = Img_origin_before_crop(Highbegin:cropHeight, Widthbegin:cropWidth);
        Img1_1 = bwareaopen(Img_origin,500,8);
        Img2_1 = imcomplement(Img1_1);
        Img3_1 = bwareaopen(Img2_1,100,8);
        Img4_1 = imcomplement(Img3_1);
        [centers, radii, metric] = imfindcircles(Img4_1,[20 45],...
            'ObjectPolarity','dark','Sensitivity',0.9);
        centers_ascend = sortrows(centers,1);
        combination    = [centers, radii];
        radii_ascend   = sortrows(combination,1);
        radii_left_all(i)  = radii_ascend(1,3);
        radii_right_all(i) = radii_ascend(2,3);
        if i==number(k)
            radii_left_average  = mean(radii_left_all);
            radii_right_average = mean(radii_right_all);
            bigradius    = max(radii_left_average, radii_right_average);
            dradius      = abs(radii_left_average - radii_right_average);
            errorradius(k) = dradius/bigradius;
        end
        imge_diameter       = sum(radii);
        real_diameter       = imge_diameter*factor;
        set(gca,'YDir','reverse');
        viscircles(centers_ascend, radii,'EdgeColor','b');
        centers_all{i}      = centers_ascend;
        real_diameter_all(i)= real_diameter;
        plot(centers_ascend(1,1),centers_ascend(1,2),'*'); hold on
        plot(centers_ascend(2,1),centers_ascend(2,2),'*');
        axis equal;
    end

    % ════════════════════════════════════════════════════════════════════
    % STEP 4 — post-separation: track droplets in frames after bridge breaks
    % ════════════════════════════════════════════════════════════════════
    post_tracking_ok = false;
    if ~isempty(separation_frame)
        num_post_frames = min(interval_post, lengthFiles - separation_frame);
        if num_post_frames >= 1
            post_tracking_ok = true;
            for j = 0:num_post_frames
                frame_idx = separation_frame + j;
                Img_post_raw = imread(strcat(path,dirname,'\add1\',...
                    files(frame_idx+gap).name));
                Img_post = Img_post_raw(Highbegin:cropHeight, Widthbegin:cropWidth);
                Img_p1 = bwareaopen(Img_post,500,8);
                Img_p2 = imcomplement(Img_p1);
                Img_p3 = bwareaopen(Img_p2,100,8);
                Img_p4 = imcomplement(Img_p3);
                [centers_p, radii_p, ~] = imfindcircles(Img_p4,[15 60],...
                    'ObjectPolarity','dark','Sensitivity',0.92);
                % after separation the droplets may be larger (merged mass
                % could split into unequal fragments) so search range is wider
                if size(centers_p,1) >= 2
                    comb_p   = sortrows([centers_p, radii_p], 1);
                    centers_post{j+1}       = comb_p(:,1:2);
                    radii_post_left_all(j+1)  = comb_p(1,3);
                    radii_post_right_all(j+1) = comb_p(2,3);
                    % visualise post-separation tracking in a separate figure
                    figure(100+k), imshow(Img_p4);
                    set(gca,'YDir','reverse');
                    viscircles(comb_p(:,1:2), comb_p(:,3),'EdgeColor','r');
                    hold on;
                    plot(comb_p(1,1),comb_p(1,2),'rs','MarkerSize',8);
                    plot(comb_p(2,1),comb_p(2,2),'rs','MarkerSize',8);
                    title(sprintf('Post-separation frame %d',j));
                    axis equal;
                else
                    warning('Frame %d: fewer than 2 droplets detected after separation.',...
                        frame_idx);
                    post_tracking_ok = false;
                    break
                end
            end
        end
    end

    % ════════════════════════════════════════════════════════════════════
    % STEP 5 — pre-collision velocity & We, B  (original logic unchanged)
    % ════════════════════════════════════════════════════════════════════
    for n = number(k)-1:-1:number(k)-interval
        left  = centers_all{1,number(k)}(1,:) - centers_all{1,n}(1,:);
        right = centers_all{1,number(k)}(2,:) - centers_all{1,n}(2,:);
        left_mold  = norm(left);
        right_mold = norm(right);
        left_right = dot(left,right);
        alpha = acos(left_right/(left_mold*right_mold));

        velocity_left  = factor*left_mold  / ((number(k)-n)*t);
        velocity_left_all(number(k)-n,:)  = velocity_left;
        velocity_left_average  = mean(velocity_left_all);

        velocity_right = factor*right_mold / ((number(k)-n)*t);
        velocity_right_all(number(k)-n,:) = velocity_right;
        velocity_right_average = mean(velocity_right_all);

        relative_velocity = sqrt(velocity_left^2 + velocity_right^2 ...
            - 2*velocity_left*velocity_right*cos(alpha));
        relative_velocity_all(number(k)-n,:) = relative_velocity;
        relative_velocity_average = mean(relative_velocity_all);

        We = Rho*sum(real_diameter_all)*relative_velocity^2 ...
            / (sigma*length(real_diameter_all));
        We_all(number(k)-n,:) = We;
        We_average = mean(We_all);

        % impact parameter B
        time1            = centers_all{n}(1,:) - centers_all{n}(2,:);
        time1_distance   = norm(time1);
        real_time1_distance = time1_distance*factor;
        time1_base       = time1/time1_distance;
        relative_position_base = (left-right)/norm(left-right);
        gamma = acos(dot(time1_base, relative_position_base) / ...
            (norm(relative_position_base)*norm(time1_base)));
        B = (real_time1_distance*abs(sin(gamma))) / real_diameter;
        B_all(number(k)-n,:) = B;
        B_average     = mean(B_all);
        average_diameter = mean(real_diameter_all);
    end

    % ════════════════════════════════════════════════════════════════════
    % STEP 6 — post-separation velocity (mirrors pre-collision logic)
    %
    %  Reference point : separation_frame (first frame with 2 regions)
    %  Direction       : displacement FROM separation_frame TO later frame
    %                    (opposite to pre-collision convention, since now
    %                     droplets are moving AWAY from each other)
    % ════════════════════════════════════════════════════════════════════
    velocity_post_left_average   = NaN;
    velocity_post_right_average  = NaN;
    relative_velocity_post_average = NaN;

    if post_tracking_ok && numel(centers_post) >= 2
        for m = 2:numel(centers_post)   % m=1 is separation_frame itself
            % displacement from separation frame (index 1) to frame m
            left_post  = centers_post{m}(1,:) - centers_post{1}(1,:);
            right_post = centers_post{m}(2,:) - centers_post{1}(2,:);

            left_post_mold  = norm(left_post);
            right_post_mold = norm(right_post);

            % individual speeds after separation
            vel_left_post  = factor * left_post_mold  / ((m-1)*t);
            vel_right_post = factor * right_post_mold / ((m-1)*t);

            velocity_post_left_all(m-1,:)  = vel_left_post;
            velocity_post_right_all(m-1,:) = vel_right_post;

            % angle between post-separation motion vectors
            if left_post_mold > 0 && right_post_mold > 0
                dot_post   = dot(left_post, right_post);
                cos_alpha_post = dot_post / (left_post_mold * right_post_mold);
                cos_alpha_post = max(-1, min(1, cos_alpha_post)); % clamp for safety
                alpha_post = acos(cos_alpha_post);
            else
                alpha_post = 0;
            end

            % relative speed after separation (law of cosines)
            rel_vel_post = sqrt(vel_left_post^2 + vel_right_post^2 ...
                - 2*vel_left_post*vel_right_post*cos(alpha_post));
            relative_velocity_post_all(m-1,:) = rel_vel_post;
        end

        velocity_post_left_average   = mean(velocity_post_left_all);
        velocity_post_right_average  = mean(velocity_post_right_all);
        relative_velocity_post_average = mean(relative_velocity_post_all);

        fprintf('\n--- Post-separation velocities (time%d) ---\n', k);
        fprintf('  Left droplet speed  after separation : %.4f m/s\n',...
            velocity_post_left_average);
        fprintf('  Right droplet speed after separation : %.4f m/s\n',...
            velocity_post_right_average);
        fprintf('  Relative velocity   after separation : %.4f m/s\n',...
            relative_velocity_post_average);
    else
        fprintf('time%d: post-separation tracking not available.\n', k);
    end

    % ════════════════════════════════════════════════════════════════════
    % STEP 7 — store all results
    % ════════════════════════════════════════════════════════════════════

    % velocity restitution ratio  e = V_rel_post / V_rel_pre
    % (coefficient of restitution — how much relative velocity is recovered)
    if ~isnan(relative_velocity_post_average) && relative_velocity_average > 0
        velocity_ratio = relative_velocity_post_average / relative_velocity_average;
    else
        velocity_ratio = NaN;
    end

    fprintf('  V_rel_post / V_rel_pre (restitution ratio) : %.4f\n', velocity_ratio);

    name_all{k}                         = dirname;
    We_final(k,:)                       = We_average;
    B_final(k,:)                        = B_average;
    D_final(k,:)                        = average_diameter;
    velocity_left_final(k,:)            = velocity_left_average;
    velocity_right_final(k,:)           = velocity_right_average;
    relative_velocity_final(k,:)        = relative_velocity_average;
    gamma_final(k,:)                    = gamma;
    velocity_post_left_final(k,:)       = velocity_post_left_average;
    velocity_post_right_final(k,:)      = velocity_post_right_average;
    relative_velocity_post_final(k,:)   = relative_velocity_post_average;
    velocity_ratio_final(k,:)           = velocity_ratio;
end

% ════════════════════════════════════════════════════════════════════════
% STEP 8 — write Excel output (columns A–M)
% Column layout:
%   A: case name      B: We         C: B (impact param)   D: diameter
%   E: V_left (pre)   F: V_right (pre)  G: V_rel (pre)    H: gamma
%   I: size error     J: V_left (post)  K: V_right (post)  L: V_rel (post)
%   M: V_rel_post / V_rel_pre  (restitution ratio)
% ════════════════════════════════════════════════════════════════════════
excelname = 'We-B-8-8-v-out';

% write header row
headers = {'Case','We','B','D(m)','V_left_pre(m/s)','V_right_pre(m/s)',...
    'V_rel_pre(m/s)','Gamma(rad)','Size_error',...
    'V_left_post(m/s)','V_right_post(m/s)','V_rel_post(m/s)',...
    'Vpost_Vpre_ratio'};
xlswrite(strcat(path,excelname), headers,                     1, 'A1');

% write data
xlswrite(strcat(path,excelname), name_all',                   1, 'A2');
xlswrite(strcat(path,excelname), We_final,                    1, 'B2');
xlswrite(strcat(path,excelname), B_final,                     1, 'C2');
xlswrite(strcat(path,excelname), D_final,                     1, 'D2');
xlswrite(strcat(path,excelname), velocity_left_final,         1, 'E2');
xlswrite(strcat(path,excelname), velocity_right_final,        1, 'F2');
xlswrite(strcat(path,excelname), relative_velocity_final,     1, 'G2');
xlswrite(strcat(path,excelname), gamma_final,                 1, 'H2');
xlswrite(strcat(path,excelname), errorradius',                1, 'I2');
xlswrite(strcat(path,excelname), velocity_post_left_final,    1, 'J2');
xlswrite(strcat(path,excelname), velocity_post_right_final,   1, 'K2');
xlswrite(strcat(path,excelname), relative_velocity_post_final,1, 'L2');
xlswrite(strcat(path,excelname), velocity_ratio_final,        1, 'M2');
