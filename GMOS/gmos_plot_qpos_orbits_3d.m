function gmos_plot_qpos_orbits_3d(traj, style, t_i)
%GMOS_PLOT_QPOS_ORBITS_3D Plot 3D trajectories with selectable styling.
%
% Allowed calls:
%   gmos_plot_qpos_orbits_3d(traj)
%   gmos_plot_qpos_orbits_3d(traj, style)
%   gmos_plot_qpos_orbits_3d(traj, t_i)
%   gmos_plot_qpos_orbits_3d(traj, style, t_i)

    % -------------------------
    % Flexible input parsing
    % -------------------------
    if nargin < 2
        style = struct();
        t_i = [];
    elseif nargin < 3
        if isstruct(style) || isempty(style)
            t_i = [];
        else
            t_i = style;
            style = struct();
        end
    end

    if isempty(style)
        style = struct();
    end
    if ~isstruct(style)
        error('Second input must be a style struct or a scalar time t_i.');
    end

    % -------------------------
    % Defaults
    % -------------------------
    if ~isfield(style, 'mode') || isempty(style.mode)
        style.mode = 'index';
    end
    if ~isfield(style, 'lineWidth') || isempty(style.lineWidth)
        style.lineWidth = 1.0;
    end
    if ~isfield(style, 'alpha') || isempty(style.alpha)
        style.alpha = 0.25;
    end

    mode = lower(string(style.mode));
    [~, ~, M] = size(traj.Y);

    if mode == "index"
        cols = lines(max(M,1));
    elseif mode == "gray"
        if ~isfield(style, 'gray') || isempty(style.gray)
            style.gray = 0.25;
        end
        cols = repmat([style.gray style.gray style.gray], M, 1);
    elseif mode == "rgb"
        if ~isfield(style, 'rgb') || isempty(style.rgb)
            style.rgb = [0 0 0];
        end
        cols = repmat(style.rgb(:).', M, 1);
    else
        error("style.mode must be 'index', 'gray', or 'rgb'.");
    end

    figure;
    hold on; grid on; axis equal;
    view(3);
    xlabel('x [n.d]'); ylabel('y [n.d.]'); zlabel('z [n.d.]');
    %title('GMOS: propagated QPO trajectories');

    for m = 1:M
        r = squeeze(traj.Y(1:3,:,m));

        patch( ...
            'XData', [r(1,:) nan], ...
            'YData', [r(2,:) nan], ...
            'ZData', [r(3,:) nan], ...
            'FaceColor', 'none', ...
            'EdgeColor', cols(m,:), ...
            'EdgeAlpha', style.alpha, ...
            'LineWidth', style.lineWidth);
    end

    % Optional ring at fixed time t_i
    if ~isempty(t_i)
        if ~isscalar(t_i) || ~isnumeric(t_i) || ~isfinite(t_i)
            error('t_i must be a finite numeric scalar.');
        end

        if isfield(traj, 't')
            tvec = traj.t;
        elseif isfield(traj, 'tout')
            tvec = traj.tout;
        else
            error('traj must contain a time vector field, such as traj.t or traj.tout.');
        end

        [~, k_i] = min(abs(tvec - t_i));

        x_ring = squeeze(traj.Y(1,k_i,:));
        y_ring = squeeze(traj.Y(2,k_i,:));
        z_ring = squeeze(traj.Y(3,k_i,:));

        plot3(x_ring, y_ring, z_ring, '-', ...
                'LineWidth', 1.5);
    end
end