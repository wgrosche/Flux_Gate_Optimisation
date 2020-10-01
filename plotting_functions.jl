
function equalize_3d_aspect()
    limmin = minimum([ l[1] for l in [xlim(), ylim(), zlim()] ])
    limmax = maximum([ l[2] for l in [xlim(), ylim(), zlim()] ])
    xlim(limmin, limmax)
    ylim(limmin, limmax)
    zlim(limmin, limmax)
end

function plot_initial_final(g, vertex_positions, poi, poi_2,  g_2, vertex_positions_2; standalone = true, alpha = 1, titles = "Title")
    standalone && subplot(111, projection = "3d")

    for edge in edges(g)
        pos1 = vertex_positions[src(edge)]
        pos2 = vertex_positions[dst(edge)]
        gca().plot([pos1[1], pos2[1]], [pos1[2], pos2[2]], [pos1[3], pos2[3]],
            color = "black", alpha = alpha)
    end

    for edge2 in edges(g_2)
        pos12 = vertex_positions_2[src(edge2)]
        pos22 = vertex_positions_2[dst(edge2)]
        gca().plot([pos12[1], pos22[1]], [pos12[2], pos22[2]], [pos12[3], pos22[3]],
            color = "blue", alpha = alpha)
    end
    gca().plot(getindex.(poi, 1), getindex.(poi, 2), getindex.(poi, 3), "x", color = "red")
    gca().plot(getindex.(poi_2, 1), getindex.(poi_2, 2), getindex.(poi_2, 3), "x", color = "green")
    if standalone
        xlabel("x")
        ylabel("y")
        title(titles)
        # tight_layout()
        equalize_3d_aspect()
    end
end

function plot_final(g, vertex_positions, poi,  g_2, vertex_positions_2; standalone = true, alpha = 1, titles = "Title")
    standalone && subplot(111, projection = "3d")

    for edge in edges(g)
        pos1 = vertex_positions[src(edge)]
        pos2 = vertex_positions[dst(edge)]
        gca().plot([pos1[1], pos2[1]], [pos1[2], pos2[2]], [pos1[3], pos2[3]],
            color = "black", alpha = alpha)
    end

    for edge2 in edges(g_2)
        pos12 = vertex_positions_2[src(edge2)]
        pos22 = vertex_positions_2[dst(edge2)]
        gca().plot([pos12[1], pos22[1]], [pos12[2], pos22[2]], [pos12[3], pos22[3]],
            color = "blue", alpha = alpha)
    end
    gca().plot(getindex.(poi, 1), getindex.(poi, 2), getindex.(poi, 3), "x", color = "red")
    if standalone
        xlabel("x")
        ylabel("y")
        title(titles)
        # tight_layout()
        equalize_3d_aspect()
    end
end


function plot_num_sens(g, vertex_positions, poi,  g_2, vertex_positions_2; standalone = true, alpha = 1, titles = "Title")
    standalone && subplot(111, projection = "3d")

    for edge in edges(g)
        pos1 = vertex_positions[src(edge)]
        pos2 = vertex_positions[dst(edge)]
        gca().plot([pos1[1], pos2[1]], [pos1[2], pos2[2]], [pos1[3], pos2[3]],
            color = "black", alpha = alpha)
    end

    for edge2 in edges(g_2)
        pos12 = vertex_positions_2[src(edge2)]
        pos22 = vertex_positions_2[dst(edge2)]
        gca().plot([pos12[1], pos22[1]], [pos12[2], pos22[2]], [pos12[3], pos22[3]],
            color = "blue", alpha = alpha)
    end
    for i in 1:17
        gca().plot(getindex.(poi[i], 1), getindex.(poi[i], 2), getindex.(poi[i], 3), "x")
    end
    if standalone
        xlabel("x")
        ylabel("y")
        title(titles)
        # tight_layout()
        equalize_3d_aspect()
    end
end
