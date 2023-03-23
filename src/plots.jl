

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
                 prefix = "",
                 figdir = "Fig",
                 title_time_format = "yyyy-mm-dd",
                 )

    function myplot(x, t, cl = extrema(skipmissing(x)); kwargs...)
        subplots_adjust(hspace = 0.35)
        #@show size(x), typeof(x), typeof(lon), typeof(lat)

        if all(ismissing,x)
            x = fill(NaN,size(x))
        else
            x = nomissing(x,NaN);
        end

        if VersionNumber(PyPlot.matplotlib.__version__) <  v"3.3"
            pcolor(lon,lat,x'; cmap=cmap, kwargs...)
        else
            pcolor(lon,lat,x'; cmap=cmap, shading = "nearest", kwargs...)
        end

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

        fig.suptitle("Date: $(Dates.format(time[n],title_time_format))", fontsize=16)
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

        figname = joinpath(figdir,prefix * replace(basename(fname_rec),".nc" => "_" * Dates.format(time[n],"yyyy-mm-ddTHHMM") * ".png"))
        @debug figname
        savefig(figname,dpi=300)

        GC.gc()
        pyGC = pyimport("gc")
        pyGC.collect()
    end

    close(ds_cv)
    close(ds_orig)
    close(ds_rec)

    return nothing
end



function plotvec(cases,fnameavg,varnames;
                 figdir = joinpath(dirname(fnameavg),"Fig"),
                 title_time_format = "yyyy-mm-dd HH:MM",
                 prefix = "",
                 maskname = nothing,
                 cl = (-0.7, 0.7),
                 ireduce = 5,
                 which_plot = :all,
                 )

    function plotrad(hfsite,r,cl)
        pcolormesh(lon,lat,r');
        clim(cl)
        colorbar()
        plot(hfsite.lon,hfsite.lat,"mo")
        DINCAE_utils.set_aspect_ratio()
        #DINCAE_utils.plotmap()
    end

    if ireduce isa Number
        jreduce = ireduce
    else
        (ireduce,jreduce) = ireduce
    end

    mkpath(figdir)

    ds = NCDataset(fnameavg)
    lon = ds["lon"][:]
    lat = ds["lat"][:]

    fnames = [case.fname_orig for case in cases]
    fnames_cv = [case.fname_cv for case in cases]

    ds_obs = NCDataset.(fnames)
    ds_obs_cv = NCDataset.(fnames_cv)
    time = ds_obs[1]["time"][:]


    fig = figure()

    #clf(); pcolormesh(lon,lat,nomissing(ds["u"][:,:,1],NaN)')

    if maskname !== nothing
        mask = NCDataset(maskname,"r") do dsmask
            ds["mask"][:,:,:] .== 1
        end
    else
        mask = trues(length(lon),length(lat))
    end

    time_index =
        if which_plot == :all
            1:length(time)
        else
            sort(unique(reduce(vcat,listcvimages.(cases))))
        end

    #for n = 1:32
    #for n = 1:length(time)
    for n = time_index
        clf()
        println("plot ",time[n])

        u = ds["u"][:,:,n]
        v = ds["v"][:,:,n]

        u_error = ds["u_error"][:,:,n]
        v_error = ds["v_error"][:,:,n]

        P11 = u_error.^2
        P22 = v_error.^2
        P12 = ds["u_v_covar"][:,:,n]

        fig.suptitle("Date: $(Dates.format(time[n],title_time_format))", fontsize=16)

        u2 = copy(u)
        v2 = copy(v)
        u2[.!mask] .= NaN
        v2[.!mask] .= NaN


        t_error = sqrt.(u_error.^2 + v_error.^2)
        uv_norm = sqrt.(u.^2 + v.^2)
        #u[t_error .>= 0.5] .= NaN
        #v[t_error .>= 0.5] .= NaN
        #u[uv_norm .>= 1.5] .= NaN
        #v[uv_norm .>= 1.5] .= NaN

        isubplot = 0

        for k = 1:length(cases)
            #global isubplot

            hfsite = cases[k]
            varname = cases[k].varname
            #ur = extract_radvel(lon,lat,u,v,hfsite.lon,hfsite.lat,hfsite.positive_toward)

            direction_obs = DINCAE_utils.dirobs(lon,lat,hfsite.lon,hfsite.lat,hfsite.positive_toward)
            ur,sigma_ur = DINCAE.vector2_projection((u,v),(P11,P12,P22),direction_obs)


            ur_obs = nomissing(ds_obs[k][varname][:,:,n],NaN)
            ur_obs_cv = nomissing(ds_obs_cv[k][varname][:,:,n],NaN)

            ur[isnan.(ur_obs)] .= NaN
            sigma_ur[isnan.(ur_obs)] .= NaN

            #cl = extrema(filter(isfinite,ur_obs))
            isubplot += 1; subplot(length(cases),4,isubplot);
            plotrad(hfsite,ur_obs,cl);

            isubplot += 1; subplot(length(cases),4,isubplot);
            plotrad(hfsite,ur_obs_cv,cl);

            isubplot += 1; subplot(length(cases),4,isubplot);
            plotrad(hfsite,ur,cl);

            ir = 1:ireduce:length(lon)
            jr = 1:jreduce:length(lat)
            quiver(lon[ir],lat[jr],u2[ir,jr]',v2[ir,jr]')

            isubplot += 1; subplot(length(cases),4,isubplot);
            plotrad(hfsite,sigma_ur,(0,0.02));
        end

        GC.gc()
        pyGC = pyimport("gc")
        pyGC.collect()

        figname = joinpath(figdir,prefix * replace(basename(fnameavg),".nc" => "_" * Dates.format(time[n],"yyyy-mm-ddTHHMM") * ".png"))
        @debug figname
        savefig(figname,dpi=300)
    end
    close.(ds_obs)
    close(ds)
end
