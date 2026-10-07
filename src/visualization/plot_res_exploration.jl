

####################### res exploration ############################
ft_labels(x) = x 
ft_ticks(x) = x -1
ft_legends(x) = x - 1
ft_insidelabels(x) = x -2
ft_title(x) = x + 1


function plot_res_exploration(final_path::String; horiz_lines::Vector{FT} = [0.0], vert_lines::Vector{FT} = [1.0])
    data = load(string(final_path,".jld"))
    println(keys(data))
    ft = 10
    rasterization = true
    #rc("font", family="LatinModernMath", size=ft)
    rc("mathtext", fontset="stix")

    O0vec = data["input_args"]["OP0"]
    toplot = data["OPnames"]

    fig = figure(figsize=(20, 8.5))

    gs_tables = gspec.GridSpec(length(O0vec), 1)
    gs_tables[:update](bottom=0.0, top=1.0, left=0.0, right=0.06, wspace=0.0)
    axs_tables = collect([plt.subplot(get(gs_tables, (jj, 0))) for jj=0:1:length(O0vec)-1])


    gs_OP = gspec.GridSpec(length(O0vec),length(toplot)+2)
    gs_OP[:update](bottom=0.0, top=1.0, left=0.1, right=1.0, wspace=0.15, hspace = 0.1)

    axs_OP = Matrix{PyObject}(undef,length(O0vec),length(toplot))
    for ii=1:length(O0vec)
        for jj=1:length(toplot)
            axs_OP[ii,jj] = plt.subplot(get(gs_OP,(ii-1,jj-1)))
        end
    end

    axs_conv = collect([plt.subplot(get(gs_OP, (jj, length(toplot)))) for jj=0:1:length(O0vec)-1])
    axs_frank = collect([plt.subplot(get(gs_OP, (jj, length(toplot)+1))) for jj=0:1:length(O0vec)-1])


    for iO=1:length(O0vec)
        #axs_tables[iO].table(Matrix(data["input_args"]["OP0"][iO]')'[:,:] , loc = "center right", cellLoc = "center", rowLabels = data["OPnames"], colWidths = fill(1/3, 1), rowLoc="center", colLoc="center")
        axs_tables[iO].set_xticklabels([]); axs_tables[iO].set_yticklabels([]);  axs_tables[iO].tick_params(axis="both",which="both",bottom=false,top=false,left=false,right=false ); axs_tables[iO].set_ylabel(string("init ",iO))
        axs_tables[iO].table(Matrix(data["input_args"]["OP0"][iO]')'[:,:] , bbox = [0.5, 0., .5, 1.], rowLabels = data["OPnames"], cellLoc = "center", colWidths = fill(1/3, 1),rowLoc="center", colLoc="center")
        for cardinal in ["top"; "left"; "right"; "bottom"]; axs_tables[iO].spines[cardinal].set_visible(false); end        
    end

    # 
    loop_params = data["loop_params"]
    xvalues = data["input_args"]["param_loop"][1][6]=="reciprocal" ? 1 ./ data["ploop_values"][loop_params[1]] : data["ploop_values"][loop_params[1]]
    yvalues = data["input_args"]["param_loop"][2][6]=="reciprocal" ? 1 ./ data["ploop_values"][loop_params[2]] : data["ploop_values"][loop_params[2]]

    xlab = data["input_args"]["param_loop"][1][6]=="reciprocal" ? string("1 / ",replace_Greek(string(data["loop_params"][1]))) :  replace_Greek(string(data["loop_params"][1]))
    ylab = data["input_args"]["param_loop"][2][6]=="reciprocal" ? string("1 / ",replace_Greek(string(data["loop_params"][2]))) :  replace_Greek(string(data["loop_params"][2]))

    cm_pos = get_cmap(:Blues)
    cm_neg = get_cmap(:Reds_r)
    
    eps = 1e-15
    vals = zeros(length(xvalues), length(yvalues))
    @showprogress for (ik,k) in enumerate(toplot)
        for iO = 1:length(O0vec)
            axs_OP[iO,ik].set_xlim(extrema(xvalues))
            axs_OP[iO,ik].set_ylim(extrema(yvalues))
            axs_OP[iO,1].set_ylabel(ylab, fontsize = ft_labels(ft))
            axs_OP[size(axs_OP,1),ik].set_xlabel(xlab, fontsize = ft_labels(ft))

            if !all(isnan.(data["res"][k][:,:,iO]))
                #negative values first
                #vals .= min.(0.0,data["res"][k][:,:,iO]); vmin, vmax = minimum( x->isnan(x) ? Inf : x,vals), maximum(x->isnan(x) ? -Inf : x,vals); replace!(vals, 0.0=>NaN)
                vals .= min.(0.0,data["res"][k][:,:,iO]); vmin, vmax = [-1.0; 0.0]; replace!(vals, 0.0=>NaN)
                axs_OP[iO,ik].pcolormesh(xvalues,yvalues, vals' , cmap = cm_neg, rasterized=rasterization)
                cbar_neg = fig.colorbar(cms.ScalarMappable( norm=pltcolors.Normalize(vmin=vmin, vmax = vmax , clip=false), cmap=cm_neg),ax=axs_OP[iO,ik], ticks = [vmin; vmax], aspect = 7, fraction = 0.04, pad = 0.02); #anchor = (-0.8,0.0))

                # plot positive values
                #vals .= max.(-eps,data["res"][k][:,:,iO]); vmin, vmax = minimum( x->isnan(x) ? Inf : x,vals), maximum(x->isnan(x) ? -Inf : x,vals); replace!(vals, -eps=>NaN)
                vals .= max.(-eps,data["res"][k][:,:,iO]); vmin, vmax = [0.0; 1.0]; replace!(vals, -eps=>NaN)
                axs_OP[iO,ik].pcolormesh(xvalues,yvalues, vals' , cmap = cm_pos,rasterized=rasterization)
                cbar_pos = fig.colorbar(cms.ScalarMappable( norm=pltcolors.Normalize(vmin=vmin, vmax = vmax , clip=false), cmap=cm_pos),ax=axs_OP[iO,ik],  ticks = [vmin; vmax], aspect = 7, fraction = 0.04,  pad = 0.015) #anchor = (3.,0.0))
                cbar_pos.ax.tick_params(labelsize = ft_insidelabels(ft) -3 ,left=true,right=false, labelleft=true, labelright=false) 
                cbar_neg.ax.tick_params(labelsize = ft_insidelabels(ft) -3 ,left=false,right=true) 
            end

            data["input_args"]["param_loop"][1][5]=="exp10" && axs_OP[iO,ik].set_xscale("log")
            data["input_args"]["param_loop"][2][5]=="exp10" && axs_OP[iO,ik].set_yscale("log")

            #xlims = axs_OP[iO,ik].get_xlim(); for hv in horiz_lines; axs_OP[iO,ik].plot(xlims, fill(hv,2), color ="black", lw = 1.5); end; axs_OP[iO,ik].set_xlim(xlims) # position horizontal lines
            #ylims = axs_OP[iO,ik].get_ylim(); for vv in vert_lines; axs_OP[iO,ik].vlines(vv, ymin = ylims[1], ymax = ylims[2], color ="black", lw = 1.5);  end; axs_OP[iO,ik].set_ylim(ylims); # position vertical lines

            iO != length(O0vec) && axs_OP[iO,ik].set_xticklabels([])
        end
        axs_OP[1,ik].set_title(string(k),fontsize = ft_title(ft))
    end

    cm_conv = pltcolors.ListedColormap(["red"; "grey";"limegreen"])
    cm0 = get_cmap(:tab10)
    cm_rank = pltcolors.ListedColormap(collect([cm0(i)[1:3] for i=0:length(O0vec)-1]))


   

    # conv
    for iO = 1:length(O0vec)
        vals .= Int.(data["res"]["conv"][:,:,iO])
        @assert sum(isnan.(vals))==0
        @assert length(unique(vals))<=3
        axs_conv[iO].set_xlim(extrema(xvalues));axs_conv[iO].set_ylim(extrema(yvalues))
            
        axs_conv[iO].pcolormesh(xvalues, yvalues, add_rgb_to_discrete_pcolormesh(vals', cm_conv, [-1;0;1]),rasterized=rasterization)
        iO != length(O0vec) && axs_conv[iO].set_xticklabels([])
        println("init -> ",O0vec[iO],"  convergence countmap -> ", countmap(vals))

        data["input_args"]["param_loop"][1][5]=="exp10" && axs_conv[iO].set_xscale("log")
        data["input_args"]["param_loop"][2][5]=="exp10" && axs_conv[iO].set_yscale("log")


        #xlims = axs_conv[iO].get_xlim(); for hv in horiz_lines; axs_conv[iO].plot(xlims, fill(hv,2), color ="black", lw = 1.5); end; axs_conv[iO].set_xlim(xlims) # position horizontal lines
        #ylims = axs_conv[iO].get_ylim(); for vv in vert_lines; axs_conv[iO].vlines(vv, ymin = ylims[1], ymax = ylims[2], color ="black", lw = 1.5); end; axs_conv[iO].set_ylim(ylims) # position vertical lines
    end

    #old command to set norm with boundary norm -> normm = pltcolors.BoundaryNorm([-1.5; -.5; .5; 1.5], cm_conv.N+1)
    ax_cbar_conv = axs_conv[1].inset_axes([0.0,1.0,1.0,0.2])
    cbar_conv = fig.colorbar(cms.ScalarMappable(cmap=cm_conv),cax=ax_cbar_conv, orientation="horizontal",location="top", ticks = [], pad = 0.0)
    for (iv, xv) in enumerate([-1; 0;1])
        cbar_conv.ax.text(iv / 3 - 1/6 , 0.5, xv, ha="center", va="center", fontsize = ft-2)
    end
    cbar_conv.ax.tick_params( direction="in",bottom=false,top = false,left=false,right=false,labelbottom =false,labeltop = false, labelleft=false, labelright=true)
    axs_conv[length(O0vec)].set_xlabel(xlab, fontsize = ft_labels(ft))
    axs_conv[1].set_title("conv",fontsize = ft_title(ft))


    
    # frank
    for iO = 1:length(O0vec)
        vals .= data["res"]["f_ranking"][:,:,iO]
        axs_frank[iO].set_xlim(extrema(xvalues));axs_frank[iO].set_ylim(extrema(yvalues))
        @assert sum(isnan.(vals))==0
        axs_frank[iO].pcolormesh(xvalues, yvalues, add_rgb_to_discrete_pcolormesh(vals', cm_rank, collect(1:length(O0vec))), cmap = cm_rank,rasterized=rasterization)
        iO != length(O0vec) && axs_frank[iO].set_xticklabels([])

        data["input_args"]["param_loop"][1][5]=="exp10" && axs_frank[iO].set_xscale("log")
        data["input_args"]["param_loop"][2][5]=="exp10" && axs_frank[iO].set_yscale("log")

        #xlims = axs_frank[iO].get_xlim(); for hv in horiz_lines; axs_frank[iO].plot(xlims, fill(hv,2), color ="black", lw = 1.5); end; axs_frank[iO].set_xlim(xlims) # position horizontal lines
        #ylims = axs_frank[iO].get_ylim(); for vv in vert_lines; axs_frank[iO].vlines(vv, ymin = ylims[1], ymax = ylims[2], color ="black", lw = 1.5); end; axs_frank[iO].set_ylim(ylims) # position vertical lines
    end
  
    ax_cbar_rank = axs_frank[1].inset_axes([0.0,1.0,1.0,0.2])
    cbar_rank = fig.colorbar(cms.ScalarMappable(cmap=cm_rank),cax=ax_cbar_rank, orientation="horizontal",location="top", ticks = [], pad = 0.0)
    for (iv, xv) in enumerate(collect(1:length(O0vec)))
        cbar_rank.ax.text(xv / length(O0vec) .- 1/(2*length(O0vec)), 0.5, xv, ha="center", va="center", fontsize = ft-2)
    end
    cbar_rank.ax.tick_params( direction="in",bottom=false,top = false,left=false,right=false,labelbottom =false,labeltop = false, labelleft=false, labelright=true)
    axs_frank[length(O0vec)].set_xlabel(xlab, fontsize = ft_labels(ft))
    axs_frank[1].set_title("f ranking",fontsize = ft_title(ft))
    



  

    #axs_OP[1,size(axs_OP,2)].table(vcat(data["input_args"]["OP0"]'...)[:,:], loc = "top", cellLoc = "center", rowLabels = vcat([string("Init ",i) for i=1:length(O0vec)]), colLabels = data["OPnames"])
    
    fig.suptitle(final_path,fontsize = ft_title(ft)+1, x = 0.1, y = 1.08, horizontalalignment = "left" )
    gcf()
    println("plot done, now saving")
    PyPlot.savefig(string(final_path,".pdf"),bbox_inches="tight", dpi=300)
    println("saving completed \n")
    return fig

end

function add_rgb_to_discrete_pcolormesh(in, cm_conv, values::Vector{Int})
    dx,dy = size(in)
    out = zeros(dx,dy,3)


    rgb_list = [cm_conv(i)[1:3] for i=0:length(values)-1]
    for (ii,vv) in enumerate(values)
            
        for CI in findall(x->x==vv, in)
            out[CI,:] .= rgb_list[ii]
        end

    end
    
    return out
    
end
