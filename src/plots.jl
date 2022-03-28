

function PyObject(a::Array{Union{Missing},N}) where {N}
  numpy_ma = pyimport("numpy").ma
  pycall(numpy_ma.array, Any, coalesce.(a,1.), mask=ismissing.(a))
end


"""
Fixes the aspect ratio of a plot.
"""
function set_aspect_ratio()
    ax = gca()
    as = cosd(mean(ylim()))
    ax.set_aspect(1/as)
end

function plotmap(; patchcolor = [.8,.8,.8], resolution='f', grid=5)
    lon,lat,data = GeoDatasets.landseamask(resolution=resolution,grid=grid)

#    function plot()
        xl = xlim()
        yl = ylim()
        #    @show size(data)
        contourf(lon,lat,data',levels = [0.5,2],colors = [patchcolor]);
        xlim(xl)
        ylim(yl)
#    end
#    return plot
end

function plotres(case, fname_rec; transfun = (identity,identity), clim = case.clim, which_plot = :all, cb_orientation = "vertical",cmap="viridis",
                 clim_quantile = (0.,1.),
                 figdir = "Fig"
                 )

    function myplot(x, t, cl = extrema(skipmissing(x)); kwargs...)
        subplots_adjust(hspace = 0.35)
        #@show size(x), typeof(x), typeof(lon), typeof(lat)

        if all(ismissing,x)
            x = fill(NaN,size(x))
        else
            x = nomissing(x,NaN);
        end
        #pcolor(lon,lat,x'; cmap=cmap, shading = "nearest", kwargs...)
        pcolor(lon,lat,x'; cmap=cmap, kwargs...)
        set_aspect_ratio()
        title(t)

        if cl != nothing
            PyPlot.clim(cl)
        end
        colorbar(orientation = cb_orientation)
        plotmap()
        xlim(extrema(lon))
        ylim(extrema(lat))
    end

    #_plotmap = plotmap()

    PyPlot.ioff()
    mkpath(figdir)
    fname_orig = case.fname_orig
    fname_cv = case.fname_cv
    varname = case.varname

    cl = clim
    ds_cv = Dataset(fname_cv);
    ds_orig = Dataset(fname_orig);
    ds_rec = Dataset(fname_rec);
    lon = ds_orig["lon"][:]
    lat = ds_orig["lat"][:]
    time = ds_orig["time"][:]
    #fn = log10
    fn = transfun[1]

    mask = ds_cv["mask"][:,:] .== 1

    time_index =
        if which_plot == :all
            1:length(time)
        else
            listcvimages(case)
        end

    fig = figure(figsize=(10,7))

    for n = time_index
    #for n = time_index[1:1]
    #for n = 1:1
        clf()
        println("plot ",time[n])

        data_cv = ds_cv[varname][:,:,n]
        data_orig = ds_orig[varname][:,:,n]
        #data_rec = ds_rec["mean_rec"][:,:,n]
        #data_rec = ds_rec["batch_m_rec"][:,:,n]
        data_rec = ds_rec[varname][:,:,n]

        data_orig[.!mask] .= missing

        cl =
            if (clim == nothing)
                tmp = skipmissing(vcat(fn.(data_orig[mask]),fn.(data_rec[mask])))
                #extrema(tmp)
                quantile(tmp,clim_quantile)
            else
                clim
            end

        @debug begin
            @show cl,time[n]
        end

        if transfun[1] == log
            norm = matplotlib.colors.LogNorm(vmin=cl[1], vmax=cl[1])
        else
            norm = matplotlib.colors.Normalize(vmin=cl[1], vmax=cl[1])
        end

        fig.suptitle("Date: $(Dates.format(time[n],"yyyy-mm-dd"))", fontsize=16)
        hspace = 0.2
        subplot(2,2,1)
        myplot(fn.(data_orig),"(a) Original data",cl,norm=norm)

        ncv = sum(ismissing.(data_cv) .&  .!ismissing.(data_orig))
        if ncv !== 0
            subplot(2,2,2)
            myplot(fn.(data_cv),"(b) With added clouds ($ncv)", cl, norm=norm)
        end

        subplot(2,2,3)
        myplot(fn.(data_rec),"(c) DINCAE reconstruction", cl, norm=norm)

        subplot(2,2,4)
        #myplot(ds_rec["sigma_rec"][:,:,n],"σ")
        #myplot(ds_rec["batch_sigma_rec"][:,:,n],"σ")
        myplot(ds_rec[varname * "_error"][:,:,n],"(d) Exp. error std. dev.")

        figname = joinpath(figdir,replace(basename(fname_rec),".nc" => "_" * Dates.format(time[n],"yyyy-mm-dd") * ".png"))
        @debug figname
        savefig(figname,dpi=300)
    end

    close(ds_cv)
    close(ds_orig)
    close(ds_rec)

    return nothing
end
