function generate_simulation_video(softRobot, dof_with_time, time_arr, sim_params, environment, imc)
    % --- 1. 自动检测有效步数 ---
    % 检查 dof_with_time 第一行（时间戳）之后的数据，寻找最后一个非零时刻
    actual_steps = find(sum(abs(dof_with_time(2:end, :)), 1) > 0, 1, 'last');
    
    if isempty(actual_steps)
        warning('未检测到有效记录数据。请检查 sim_params.log_data 是否开启。');
        return;
    end
    
    fprintf('检测到已完成 %d 步仿真数据。开始生成视频...\n', actual_steps);

    % --- 2. 视频对象配置 ---
    video_name = sprintf('deicing_video_%s.mp4', datestr(now, 'yyyymmdd_HHMM'));
    v = VideoWriter(video_name, 'MPEG-4');
    v.FrameRate = 30; 
    v.Quality = 100;
    open(v);

    % --- 3. 初始化绘图窗口 ---
    h_fig = figure('Name', 'Auto-Detection Video Export', 'Color', 'w', 'Position', [100 100 800 600]);
    h_ax = axes('Parent', h_fig);

    % --- 4. 状态预处理 ---
    % 重新初始化参考系以保证回放连续性
    softRobot = computeSpaceParallel(softRobot);

    % 自动设置抽帧步长：如果数据量太大，自动增大步长以加快生成速度
    plot_skip = max(1, floor(actual_steps / 300)); % 目标生成约 300 帧左右的视频

    for t = 1 : plot_skip : actual_steps
        % 确保绘图在正确的窗口
        set(0, 'CurrentFigure', h_fig);
        set(h_fig, 'CurrentAxes', h_ax);
        
        % 同步 q 与 q0 消除 Bishop 框架闪烁
        softRobot.q = dof_with_time(2:end, t); 
        softRobot.q0 = softRobot.q; 
        
        cla(h_ax); 
        hold(h_ax, 'on');
        
        % 调用核心绘图函数
        plot_MultiRod(softRobot, time_arr(t), sim_params, environment, imc);
        
        % 轴范围容错处理
        try
            if isfield(sim_params, 'plot_x') && length(sim_params.plot_x) == 2
                ax_lims = [sim_params.plot_x(1), sim_params.plot_x(2), ...
                           sim_params.plot_y(1), sim_params.plot_y(2), ...
                           sim_params.plot_z(1), sim_params.plot_z(2)];
                axis(h_ax, ax_lims);
            else
                axis(h_ax, [-0.3, 0.3, -0.3, 0.3, -0.1, 0.5]);
            end
        catch
            axis(h_ax, 'tight');
        end
        
        view(h_ax, 45, 30); 
        grid(h_ax, 'on');
        title(h_ax, sprintf('Playback Time: %.4f s', time_arr(t)));
        
        drawnow;
        frame = getframe(h_fig);
        writeVideo(v, frame);
    end

    close(v);
    fprintf('视频已保存: %s\n', fullfile(pwd, video_name));
end