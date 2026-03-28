using MathTeXEngine
struct orbit_struct
    name::String
    tspan::LinRange
    syntraj::Vector{<:SVector{6}}
    NFtraj::Vector{<:SVector{6}}
    AAtraj::Vector{<:SVector{6}}
end
function figure_settings(font_size::Symbol = :large)
    axis_fontsize, legend_fontsize, text_fontsize, title_fontsize = if font_size == :large
        28, 28, 20, 28
    elseif font_size == :medium
        24, 24, 20, 24
    else
        20, 20, 18, 20
    end

    # Custom color order (matches MATLAB corder)
    corder = [
        RGBf(0, 0, 1),
        RGBf(1, 0, 0),
        RGBf(0.8, 0.8, 0),
        RGBf(0.66, 0.66, 0.66),
        RGBf(0, 0, 0),
        RGBf(1, 0.645, 0),
        RGBf(1, 0, 1),
        RGBf(0, 0.5, 0.5),
        RGBf(0, 0, 0.543),
        RGBf(0, 0.391, 0),
        RGBf(0, 1, 1),
        RGBf(0.598, 0.195, 0.797),
    ]

    # Width / aspect ratio
    width_cm = 20
    hwratio = 0.65
    width_px = width_cm / 2.54 * 96
    height_px = hwratio * width_px

    theme = Theme(
        # ── Figure ──────────────────────────────────────────────────────────
        Figure = (backgroundcolor = :white),

        # ── Fonts ────────────────────────────────────────────────────────────
        fonts = Attributes(
            :bold => texfont(:bold),
            :bolditalic => texfont(:bolditalic),
            :italic => texfont(:italic),
            :regular => texfont(:regular),
        ),

        # ── Axes ─────────────────────────────────────────────────────────────
        Axis = (
            backgroundcolor = :white,
            spinewidth = 1,
            xgridvisible = true,
            ygridvisible = true,
            xgridcolor = (:black, 0.15),
            ygridcolor = (:black, 0.15),
            xticklabelsize = axis_fontsize,
            yticklabelsize = axis_fontsize,
            xlabelsize = axis_fontsize,
            ylabelsize = axis_fontsize,
            xlabelpadding = 20,
            ylabelpadding = 20,
            xticks = WilkinsonTicks(5),
            yticks = WilkinsonTicks(3),
            titlesize = title_fontsize,
            titlefont = :regular,
        ),

        # ── 3-D Axes ─────────────────────────────────────────────────────────
        Axis3 = (
            backgroundcolor = :white,
            spinewidth = 1,
            xgridvisible = true,
            ygridvisible = true,
            zgridvisible = true,
            xticklabelsize = axis_fontsize,
            yticklabelsize = axis_fontsize,
            zticklabelsize = axis_fontsize,
            xlabelsize = axis_fontsize,
            ylabelsize = axis_fontsize,
            zlabelsize = axis_fontsize,
            xlabeloffset = 80,
            ylabeloffset = 100,
            zlabeloffset = 120,
            xticks = WilkinsonTicks(4),
            yticks = WilkinsonTicks(4),
            zticks = WilkinsonTicks(4),
            titlesize = title_fontsize,
            titlefont = :regular,
            xspinecolor_1 = :black,
            xspinecolor_2 = :transparent,
            xspinecolor_3 = :transparent,
            yspinecolor_1 = :black,
            yspinecolor_2 = :transparent,
            yspinecolor_3 = :transparent,
            zspinecolor_1 = :black,
            zspinecolor_2 = :transparent,
            zspinecolor_3 = :transparent,
        ),

        # ── Lines ────────────────────────────────────────────────────────────
        Lines = (linewidth = 2, cycle = Cycle([:color], covary = true)),
        Stairs = (linewidth = 1.5, cycle = Cycle([:color], covary = true)),
        Scatter = (markersize = 10, cycle = Cycle([:color], covary = true)),

        # ── Legend ───────────────────────────────────────────────────────────
        Legend = (
            labelsize = legend_fontsize,
            framevisible = true,
            backgroundcolor = :white,
        ),

        # ── Colorbar ─────────────────────────────────────────────────────────
        Colorbar = (labelsize = legend_fontsize, ticklabelsize = legend_fontsize),

        # ── Text / annotations ───────────────────────────────────────────────
        Text = (fontsize = text_fontsize,),

        # ── Color palette ────────────────────────────────────────────────────
        palette = (color = corder,),
        textlabelsize = text_fontsize,
    )

    set_theme!(theme)
    return theme
end

using GLMakie, StaticArrays, LaTeXStrings

# ============================================================
# plot_orbit_phasespace
#   Left:  animated 3D synodic trajectory
#   Right: q2 vs p2  (top),  q3 vs p3  (bottom)
# ============================================================
function plot_orbit_phasespace(orb::orbit_struct,
                                transformations;
                                fps::Int     = 30,
                                mp4_path     = nothing,
                                color        = :royalblue)

    traj  = orb.syntraj
    NF    = orb.NFtraj
    N     = length(traj)
    ts    = collect(orb.tspan)

    x  = getindex.(traj, 1)
    y  = getindex.(traj, 2)
    z  = getindex.(traj, 3)

    q2 = getindex.(NF, 3)
    p2 = getindex.(NF, 4)
    q3 = getindex.(NF, 5)
    p3 = getindex.(NF, 6)

    xL1 = transformations.recentering[1]

    # ── Figure ──────────────────────────────────────────────
    set_theme!(figure_settings(:small))
    fig = Figure(size=(1600, 800), figure_padding=(80, 20, 20, 20))

    ax3d = Axis3(fig[1:2, 1];
        xlabel=L"x [N.D]", ylabel=L"y [N.D]", zlabel=L"z [N.D]",
        aspect=:data,
        title=orb.name,
        azimuth   = π/3,
        elevation = 0.15π,
    )

    ax_q2p2 = Axis(fig[1, 2];
        xlabel=L"\bar{q}_2", ylabel=L"\bar{p}_2",
        title=L"\bar{q}_2\ \mathrm{vs}\ \bar{p}_2",
    )
    ax_q3p3 = Axis(fig[2, 2];
        xlabel=L"\bar{q}_3", ylabel=L"\bar{p}_3",
        title=L"\bar{q}_3\ \mathrm{vs}\ \bar{p}_3",
    )

    colsize!(fig.layout, 1, Relative(0.6))
    colsize!(fig.layout, 2, Relative(0.4))

    # ghost traces (faded full orbit for reference)
    lines!(ax3d,    x,  y,  z;  color=(:black, 0.15), linewidth=1, transparency=true)
    lines!(ax_q2p2, q2, p2;     color=(:black, 0.15), linewidth=1, transparency=true)
    lines!(ax_q3p3, q3, p3;     color=(:black, 0.15), linewidth=1, transparency=true)

    # L1 point
    scatter!(ax3d, [xL1], [0.0], [0.0];
        color=:red, markersize=12, overdraw=true)

    ax_q2p2.aspect = DataAspect()
    ax_q3p3.aspect = DataAspect()

    # ── Observables ─────────────────────────────────────────
    tr3d   = Observable(Point3f[])
    trq2p2 = Observable(Point2f[])
    trq3p3 = Observable(Point2f[])
    dot3d  = Observable(Point3f[(x[1], y[1], z[1])])

    lines!(ax3d,    tr3d;   color=color, linewidth=2)
    lines!(ax_q2p2, trq2p2; color=color, linewidth=1.5)
    lines!(ax_q3p3, trq3p3; color=color, linewidth=1.5)
    scatter!(ax3d, dot3d;   color=:black, markersize=10)

    # ── Update ──────────────────────────────────────────────
    az_base = π/3
    function update!(i)
        tr3d[]   = Point3f.(x[1:i],  y[1:i],  z[1:i])
        trq2p2[] = Point2f.(q2[1:i], p2[1:i])
        trq3p3[] = Point2f.(q3[1:i], p3[1:i])
        dot3d[]  = [Point3f(x[i], y[i], z[i])]
        ax3d.azimuth[] = az_base + (π/4) * sin(2π * (i-1) / N)
    end

    # ── Record / display ─────────────────────────────────────
    frames = 1:N
    if !isnothing(mp4_path)
        record(fig, mp4_path, frames; framerate=fps) do i
            update!(i)
        end
        println("Saved: $mp4_path")
    end

    update!(N)
    display(fig)

    return fig
end


# ============================================================
# plot_orbit_phasespace_2d
#   Left:  animated x-y trajectory
#   Right: q2 vs p2  (top),  q3 vs p3  (bottom)
# ============================================================
function plot_orbit_phasespace_2d(orb::orbit_struct,
                                   transformations;
                                   fps::Int  = 30,
                                   mp4_path  = nothing,
                                   color     = :royalblue)

    traj = orb.syntraj
    NF   = orb.NFtraj
    N    = length(traj)

    x  = getindex.(traj, 1)
    y  = getindex.(traj, 2)

    q2 = getindex.(NF, 3)
    p2 = getindex.(NF, 4)
    q3 = getindex.(NF, 5)
    p3 = getindex.(NF, 6)

    xL1 = transformations.recentering[1]

    # ── Figure ──────────────────────────────────────────────
    set_theme!(figure_settings(:small))
    fig = Figure(size=(1600, 800), figure_padding=(40, 20, 20, 20))

    ax_xy = Axis(fig[1:2, 1];
        xlabel=L"x [N.D]", ylabel=L"y [N.D]",
        title=orb.name,
        aspect=DataAspect(),
    )
    ax_q2p2 = Axis(fig[1, 2];
        xlabel=L"\bar{q}_2", ylabel=L"\bar{p}_2",
        title=L"\bar{q}_2\ \mathrm{vs}\ \bar{p}_2",
        aspect=DataAspect(),
    )
    ax_q3p3 = Axis(fig[2, 2];
        xlabel=L"\bar{q}_3", ylabel=L"\bar{p}_3",
        title=L"\bar{q}_3\ \mathrm{vs}\ \bar{p}_3",
        aspect=DataAspect(),
    )

    colsize!(fig.layout, 1, Relative(0.6))
    colsize!(fig.layout, 2, Relative(0.4))

    # ghost traces
    lines!(ax_xy,   x,  y;   color=(:black, 0.15), linewidth=1, transparency=true)
    lines!(ax_q2p2, q2, p2;  color=(:black, 0.15), linewidth=1, transparency=true)
    lines!(ax_q3p3, q3, p3;  color=(:black, 0.15), linewidth=1, transparency=true)

    # L1 point
    scatter!(ax_xy, [xL1], [0.0]; color=:red, markersize=12, overdraw=true)

    # ── Observables ─────────────────────────────────────────
    tr_xy  = Observable(Point2f[])
    trq2p2 = Observable(Point2f[])
    trq3p3 = Observable(Point2f[])
    dot_xy = Observable(Point2f[(x[1], y[1])])

    lines!(ax_xy,   tr_xy;  color=color, linewidth=2)
    lines!(ax_q2p2, trq2p2; color=color, linewidth=1.5)
    lines!(ax_q3p3, trq3p3; color=color, linewidth=1.5)
    scatter!(ax_xy, dot_xy; color=:black, markersize=10, overdraw=true)

    # ── Update ──────────────────────────────────────────────
    function update!(i)
        tr_xy[]  = Point2f.(x[1:i],  y[1:i])
        trq2p2[] = Point2f.(q2[1:i], p2[1:i])
        trq3p3[] = Point2f.(q3[1:i], p3[1:i])
        dot_xy[] = [Point2f(x[i], y[i])]
    end

    # ── Record / display ─────────────────────────────────────
    if !isnothing(mp4_path)
        record(fig, mp4_path, 1:N; framerate=fps) do i
            update!(i)
        end
        println("Saved: $mp4_path")
    end

    update!(N)
    display(fig)

    return fig
end


# ============================================================
# plot_orbit_actionangle_2d
#   Left:  animated x-y trajectory
#   Right: I2(t), θ2(t), I3(t), θ3(t)
# ============================================================
function plot_orbit_actionangle_2d(orb::orbit_struct,
                                    transformations;
                                    fps::Int  = 30,
                                    mp4_path  = nothing,
                                    color     = :royalblue)

    traj = orb.syntraj
    AA   = orb.AAtraj
    N    = length(traj)
    ts   = collect(orb.tspan)

    x  = getindex.(traj, 1)
    y  = getindex.(traj, 2)

    I2 = getindex.(AA, 3)
    θ2 = mod.(getindex.(AA, 4), 2π)
    I3 = getindex.(AA, 5)
    θ3 = mod.(getindex.(AA, 6), 2π)

    xL1 = transformations.recentering[1]

    # ── Figure ──────────────────────────────────────────────
    set_theme!(figure_settings(:small))
    fig = Figure(size=(1600, 900), figure_padding=(40, 20, 20, 20))

    ax_xy = Axis(fig[1:4, 1];
        xlabel=L"x [N.D]", ylabel=L"y [N.D]",
        title=orb.name,
        aspect=DataAspect(),
    )

    function side_ax(row, title_str)
        Axis(fig[row, 2];
            title=title_str,
            xlabel=L"t\ \mathrm{[N.D]}",
        )
    end

    ax_I2 = side_ax(1, L"\hat{I}_2(t)")
    ax_θ2 = side_ax(2, L"\theta_2(t)\ \mathrm{[rad]}")
    ax_I3 = side_ax(3, L"\hat{I}_3(t)")
    ax_θ3 = side_ax(4, L"\theta_3(t)\ \mathrm{[rad]}")

    colsize!(fig.layout, 1, Relative(0.6))
    colsize!(fig.layout, 2, Relative(0.4))

    # ghost traces
    lines!(ax_xy,  x,  y;   color=(:black, 0.15), linewidth=1, transparency=true)
    lines!(ax_I2, ts, I2;   color=(:black, 0.15), linewidth=1, transparency=true)
    lines!(ax_θ2, ts, θ2;   color=(:black, 0.15), linewidth=1, transparency=true)
    lines!(ax_I3, ts, I3;   color=(:black, 0.15), linewidth=1, transparency=true)
    lines!(ax_θ3, ts, θ3;   color=(:black, 0.15), linewidth=1, transparency=true)

    # L1 point
    scatter!(ax_xy, [xL1], [0.0]; color=:red, markersize=12, overdraw=true)

    # ── Observables ─────────────────────────────────────────
    tr_xy  = Observable(Point2f[])
    trI2   = Observable(Point2f[])
    trθ2   = Observable(Point2f[])
    trI3   = Observable(Point2f[])
    trθ3   = Observable(Point2f[])
    dot_xy = Observable(Point2f[(x[1], y[1])])

    lines!(ax_xy,  tr_xy;  color=color, linewidth=2)
    lines!(ax_I2,  trI2;   color=color, linewidth=1.5)
    lines!(ax_θ2,  trθ2;   color=color, linewidth=1.5)
    lines!(ax_I3,  trI3;   color=color, linewidth=1.5)
    lines!(ax_θ3,  trθ3;   color=color, linewidth=1.5)
    scatter!(ax_xy, dot_xy; color=:black, markersize=10, overdraw=true)

    # ── Update ──────────────────────────────────────────────
    function update!(i)
        tr_xy[]  = Point2f.(x[1:i],  y[1:i])
        trI2[]   = Point2f.(ts[1:i], I2[1:i])
        trθ2[]   = Point2f.(ts[1:i], θ2[1:i])
        trI3[]   = Point2f.(ts[1:i], I3[1:i])
        trθ3[]   = Point2f.(ts[1:i], θ3[1:i])
        dot_xy[] = [Point2f(x[i], y[i])]
    end

    # ── Record / display ─────────────────────────────────────
    if !isnothing(mp4_path)
        record(fig, mp4_path, 1:N; framerate=fps) do i
            update!(i)
        end
        println("Saved: $mp4_path")
    end

    update!(N)
    display(fig)

    return fig
end
# ============================================================
# plot_orbit_actionangle
#   Left:  animated 3D synodic trajectory
#   Right: I2(t), θ2(t), I3(t), θ3(t)
# ============================================================
function plot_orbit_actionangle(orb::orbit_struct,
                                 transformations;
                                 fps::Int    = 30,
                                 mp4_path    = nothing,
                                 color       = :royalblue)

    traj = orb.syntraj
    AA   = orb.AAtraj
    N    = length(traj)
    ts   = collect(orb.tspan)

    x  = getindex.(traj, 1)
    y  = getindex.(traj, 2)
    z  = getindex.(traj, 3)

    I2 = getindex.(AA, 3)
    θ2 = mod.(getindex.(AA, 4), 2π)
    I3 = getindex.(AA, 5)
    θ3 = mod.(getindex.(AA, 6), 2π)

    xL1 = transformations.recentering[1]

    # ── Figure ──────────────────────────────────────────────
    set_theme!(figure_settings(:small))
    fig = Figure(size=(1600, 900), figure_padding=(80, 20, 20, 20))

    ax3d = Axis3(fig[1:4, 1];
        xlabel=L"x [N.D]", ylabel=L"y [N.D]", zlabel=L"z [N.D]",
        aspect=:data,
        title=orb.name,
        azimuth   = π/3,
        elevation = 0.15π,
    )

    function side_ax(row, title_str)
        Axis(fig[row, 2];
            title=title_str,
            xlabel=L"t\ \mathrm{[N.D]}",
        )
    end

    ax_I2 = side_ax(1, L"\hat{I}_2(t)")
    ax_θ2 = side_ax(2, L"\theta_2(t)\ \mathrm{[rad]}")
    ax_I3 = side_ax(3, L"\hat{I}_3(t)")
    ax_θ3 = side_ax(4, L"\theta_3(t)\ \mathrm{[rad]}")

    colsize!(fig.layout, 1, Relative(0.6))
    colsize!(fig.layout, 2, Relative(0.4))

    # ghost traces
    lines!(ax3d,  x,   y,   z;  color=(:black, 0.15), linewidth=1, transparency=true)
    lines!(ax_I2, ts, I2;       color=(:black, 0.15), linewidth=1, transparency=true)
    lines!(ax_θ2, ts, θ2;       color=(:black, 0.15), linewidth=1, transparency=true)
    lines!(ax_I3, ts, I3;       color=(:black, 0.15), linewidth=1, transparency=true)
    lines!(ax_θ3, ts, θ3;       color=(:black, 0.15), linewidth=1, transparency=true)

    # L1 point
    scatter!(ax3d, [xL1], [0.0], [0.0];
        color=:red, markersize=12, overdraw=true)

    # ── Observables ─────────────────────────────────────────
    tr3d  = Observable(Point3f[])
    trI2  = Observable(Point2f[])
    trθ2  = Observable(Point2f[])
    trI3  = Observable(Point2f[])
    trθ3  = Observable(Point2f[])
    dot3d = Observable(Point3f[(x[1], y[1], z[1])])

    lines!(ax3d,  tr3d;  color=color, linewidth=2)
    lines!(ax_I2, trI2;  color=color, linewidth=1.5)
    lines!(ax_θ2, trθ2;  color=color, linewidth=1.5)
    lines!(ax_I3, trI3;  color=color, linewidth=1.5)
    lines!(ax_θ3, trθ3;  color=color, linewidth=1.5)
    scatter!(ax3d, dot3d; color=:black, markersize=10)

    # ── Update ──────────────────────────────────────────────
    az_base = π/3
    function update!(i)
        tr3d[]  = Point3f.(x[1:i],   y[1:i],  z[1:i])
        trI2[]  = Point2f.(ts[1:i],  I2[1:i])
        trθ2[]  = Point2f.(ts[1:i],  θ2[1:i])
        trI3[]  = Point2f.(ts[1:i],  I3[1:i])
        trθ3[]  = Point2f.(ts[1:i],  θ3[1:i])
        dot3d[] = [Point3f(x[i], y[i], z[i])]
        ax3d.azimuth[] = az_base + (π/4) * sin(2π * (i-1) / N)
    end

    # ── Record / display ─────────────────────────────────────
    frames = 1:N
    if !isnothing(mp4_path)
        record(fig, mp4_path, frames; framerate=fps) do i
            update!(i)
        end
        println("Saved: $mp4_path")
    end

    update!(N)
    display(fig)

    return fig
end