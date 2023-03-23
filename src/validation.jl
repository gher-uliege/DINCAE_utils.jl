
"""
    rms(x,y)

Root mean square difference between `x` and `y`
"""
rms(x,y) = sqrt(mean((x - y).^2))


"""
    crms(x,y)

Centred root mean square difference between `x` and `y`
"""
crms(x,y) = rms(x .- mean(x),y .- mean(y))

bias(x,y) = mean(x) - mean(y)


function loadbatch_scalar(case,fname)
    ds = Dataset(fname);
    lon = ds["lon"][:]
    lat = ds["lat"][:]

    # all data
    ntest = ds.dim["time"]

    if haskey(case,:ntest)
        if get(case,:ntest,"") == "last-50"
            ntest = 50
        end
    end

    @show ntest

    mean_varname = case.varname
    sigma_varname = case.varname * "_error"

    batch_m_rec = ds[mean_varname][:,:,end-ntest+1:end]

    batch_sigma_rec =
        if sigma_varname != nothing
            ds[sigma_varname][:,:,end-ntest+1:end]
        else
            @warn "no error estimate found"
            zeros(size(batch_m_rec))
        end


    # add mean
    if haskey(ds,"batch_m_rec")
        # old files
        meandata = ds["meandata"][:]
        batch_m_rec = batch_m_rec .+ meandata
    end

    close(ds)

    ntest = size(batch_m_rec,3)

    ds = Dataset(case.fname_cv);
    batch_m_in = ds[case.varname][:,:,end-ntest+1:end]
    mask = ds["mask"][:];
    close(ds)

    ds = Dataset(case.fname_orig);
    batch_m_true = ds[case.varname][:,:,end-ntest+1:end]
    close(ds)

    return lon,lat,batch_m_true,batch_m_in,batch_m_rec,batch_sigma_rec,mask
end


"""
    az = azimuth(lat1,lon1,lat2,lon2)

Compute azimuth, i.e. the angle at (`lat1`,`lon1`) between the point (`lat2`,`lon2`) and the North, counted clockwise starting from the North.
The units of all input and output parameters are degrees.

```
          North
            ↑
            | .
            |   . az
(lat1,lon1) +   .
             ╲ ↙
              ╲
               ╲
                * (lat2,lon2)
```
"""
function azimuth(lat1,lon1,lat2,lon2)
    # https://en.wikipedia.org/w/index.php?title=Azimuth&oldid=750059816#Calculating_azimuth

    Δλ = π/180 * (lon2 - lon1)
    ϕ1 = π/180 * lat1
    ϕ2 = π/180 * lat2

    α = atan(sin(Δλ), cos(ϕ1)*tan(ϕ2) - sin(ϕ1)*cos(Δλ))
    return 180/π * α
end


function dirobs(lon,lat,sitelon,sitelat,positive_toward)
    direction_obs =
        if lat isa Vector
            azimuth.(lat',lon,sitelat,sitelon)
        else
            azimuth.(lat,lon,sitelat,sitelon)
        end

    # for WERA positve is away from the HF radar site, so add 180 to direction
    if !positive_toward
        direction_obs .+= 180
    end
    return direction_obs
end

function loadradials(case,fname)
    ds = Dataset(fname);
    lon = ds["lon"][:]
    lat = ds["lat"][:]

    # all data
    n = 1:ds.dim["time"]

    @show n

    mean_varname = case.varname
    sigma_varname = case.varname * "_error"

    uname,vname = case.vectorname
    u = ds[uname][:,:,n]
    v = ds[vname][:,:,n]

    P11 = ds[uname * "_error"][:,:,n].^2
    P22 = ds[vname * "_error"][:,:,n].^2
    P12 = ds[uname * "_" * vname * "_covar"][:,:,n]

    direction_obs = dirobs(lon,lat,case.lon,case.lat,case.positive_toward)[:,:,1:1]
    batch_m_rec, batch_sigma_rec = DINCAE.vector2_projection((u,v),(P11,P12,P22),direction_obs)

    ds = Dataset(case.fname_cv);
    batch_m_in = ds[case.varname][:,:,n]
    mask = ds["mask"][:];
    close(ds)

    ds = Dataset(case.fname_orig);
    batch_m_true = ds[case.varname][:,:,n]
    close(ds)

    return lon,lat,batch_m_true,batch_m_in,batch_m_rec,batch_sigma_rec,mask
end


function loadbatch(case,fname)
    if hasproperty(case,:vectorname)
        loadradials(case,fname)
    else
        loadbatch_scalar(case,fname)
    end
end


function summary(case,fname;
                 fnamesummary = replace(fname,".nc" => "-" * case.varname * ".json"),
)
    if isfile(fnamesummary)
    #if false
        summary = JSON.parse(read(fnamesummary,String))
        return summary
    else
        lon,lat,batch_m_true,batch_m_in,batch_m_rec,batch_sigma_rec,mask = loadbatch(case,fname)
        @show size(batch_m_true)
        mm = (ismissing.(batch_m_in) .&
            .!ismissing.(batch_m_true) .&
            .!ismissing.(batch_m_rec) .&
            reshape(mask .== 1,(size(mask,1),size(mask,2),1)));
        #mm = ismissing.(batch_m_in) .& .!ismissing.(batch_m_true)
        m_true = batch_m_true[mm]
        m_rec = batch_m_rec[mm]
        sigma_rec = batch_sigma_rec[mm];
        rms = sqrt(mean((m_true - m_rec).^2))
        bias = mean(m_true - m_rec)

        diff = abs.(m_true-m_rec);
        min_diff = minimum(diff)
        max_diff = maximum(diff)
        q10_diff,q90_diff = quantile(diff,(0.1,0.9))

        summary = Dict(
            "cvrms" => rms,
            "cvbias" => bias,
            "std_true" => std(m_true),
            "std_rec" => std(m_rec),
            "cor" => cor(m_true,m_rec),
            "cvcrms" => crms(m_true,m_rec),
            "number" => sum(mm),
            "min_diff" => min_diff,
            "max_diff" => max_diff,
            "q10_diff" => q10_diff,
            "q90_diff" => q90_diff,
        )

        for i = 1:3
            summary["$i-sigma"] = mean(abs.(m_true-m_rec) .<  i*sigma_rec)
        end

        open(fnamesummary,"w") do f
            JSON.print(f, summary)
        end
        return summary
    end
end

cvrms(case,fname; kwargs...) = Float32(summary(case,fname; kwargs...)["cvrms"])

function errstat(io,case,fname::AbstractString; figprefix = replace(fname,".nc" => ""))
    println(io,fname)

    lon,lat,batch_m_true,batch_m_in,batch_m_rec,batch_sigma_rec,mask = loadbatch(case,fname)

    errstat(io,case,batch_m_true,batch_m_in,batch_m_rec,batch_sigma_rec,figprefix)
end

function errstat(io,case,batch_m_true,batch_m_in,batch_m_rec,batch_sigma_rec,figprefix)
    diff = abs.(batch_m_true-batch_m_rec);


    # all
    #mm = .!(ismissing.(diff) .| ismissing.(batch_sigma_rec));
    # only CV
    mm = ismissing.(batch_m_in) .& .!ismissing.(batch_m_true)
    mm = ismissing.(batch_m_in) .& .!ismissing.(batch_m_true) .& .!ismissing.(batch_m_rec)
    println(io,"only CV points")

    m_true = batch_m_true[mm]
    m_rec = batch_m_rec[mm];
    sigma_rec = batch_sigma_rec[mm];

    println(io,"Number of CV points ",sum(mm))
    println(io,"Number of total data points ",sum(.!ismissing.(batch_m_true)))
    println(io,"RMS ",sqrt(mean((m_true - m_rec).^2)))

    for i = 1:3
        println(io,"$i-sigma: ",mean(abs.(m_true-m_rec) .<  i*sigma_rec),"  ",2*cdf(Normal(),i)-1)
    end

    x = (m_true-m_rec) ./ sigma_rec;

    if false
        clf()
        scatter(m_true,m_rec,10,sigma_rec; cmap = "jet")
        datarange = (min(minimum(m_true),minimum(m_rec)), max(maximum(m_true),maximum(m_rec)))
        xlim(datarange)
        ylim(datarange)
        plot([datarange[1],datarange[end]],[datarange[1],datarange[end]],"k--")
        xlabel("true $(case.varname)")
        ylabel("reconstructed $(case.varname)")
        axis("equal")
        colorbar()
        savefig(figprefix * "-scatter-err.png",dpi=300)
    end

    if true
        clf();
        pp = -10:0.1:10
        PyPlot.plt.hist(x,100, density = true, label = "scaled errors")
        plot(pp,pdf.(Normal(0,1),pp), label = "Normal distribution")
        xlim(pp[1],pp[end])
        ylim(0,0.5)
        legend()
        savefig(figprefix * "-pdf-err.png",dpi=300)
        savefig(figprefix * "-pdf-err.png")
        @show fit(Normal, x)
    end

end


function loadrec(fnamerec)
    Dataset(fnamerec) do ds
        if haskey(ds,"lon")
            glon = nomissing(ds["lon"][:])
            glat = nomissing(ds["lat"][:])
            gtime = nomissing(ds["time"][:])
            gu = nomissing(ds["u"][:],NaN)
            gv = nomissing(ds["v"][:],NaN)
            return glon,glat,gtime,gu,gv
        elseif haskey(ds,"LON")
            glon = nomissing(ds["LON"][:])
            glat = nomissing(ds["LAT"][:])
            gtime = nomissing(ds["time"][:])
            gu = nomissing(ds["U"][:],NaN)
            gv = nomissing(ds["V"][:],NaN)
            return glon,glat,gtime,gu,gv
        else

            glon = nomissing(ds["LONGITUDE"][:])
            glat = nomissing(ds["LATITUDE"][:])
            gtime = nomissing(ds["TIME"][:]) # .+ Dates.Hour(6)
            gu = nomissing(ds["EWCT"][:],NaN)
            gv = nomissing(ds["NSCT"][:],NaN)
            return glon,glat,gtime,gu,gv
        end

    end
end

function nanrms(a,b)
    d = a - b
    m = isfinite.(d)
    sqrt(sum(d[m].^2)/sum(m))
end

const time0 = DateTime(2000,1,1)
tonum(dt) = Dates.value(dt - time0)

function validate_drifter(
    fnamerec::AbstractString, (lon,lat,time),(u,v),ids;
    plotting = false, refcase = false,
    label = nothing,
    used_ids = sort(unique(ids)),
    figdir = joinpath(dirname(fnamerec),"Fig"),
    drifter_names =  nothing,
    exclude = nothing)


    ndrifter = length(used_ids)
    RMSstat_drifter = zeros(ndrifter)
    RMSstat_drifter_u = zeros(ndrifter)
    RMSstat_drifter_v = zeros(ndrifter)
    len_drifter = zeros(ndrifter)

    glon,glat,gtime,gu,gv = loadrec(fnamerec)

    itpu = extrapolate(interpolate((glon,glat,tonum.(gtime)),gu,Gridded(Linear())),NaN);
    itpv = extrapolate(interpolate((glon,glat,tonum.(gtime)),gv,Gridded(Linear())),NaN);

    ui = itpu.(lon,lat,tonum.(time));
    vi = itpv.(lon,lat,tonum.(time));

    if exclude != nothing
        ui[exclude] .= NaN
        vi[exclude] .= NaN
    end

    for idrifter = 1:length(used_ids)
        sel = ids .== used_ids[idrifter]

        RMSstat_drifter[idrifter] = sqrt(nanrms(ui[sel],u[sel])^2 +  nanrms(vi[sel],v[sel])^2)
        RMSstat_drifter_u[idrifter] = nanrms(ui[sel],u[sel])
        RMSstat_drifter_v[idrifter] = nanrms(vi[sel],v[sel])
        len_drifter[idrifter]  = sum(isfinite.(ui[sel]) .& isfinite.(u[sel]) .& isfinite.(vi[sel]) .& isfinite.(v[sel]))
        #@show len_drifter[idrifter]

        tmp = u[sel]

        if plotting
            fig=figure(idrifter)
            subplot(2,1,1)
#            if refcase
                subplot(2,1,1)
                plot(time[sel],u[sel]; label = "obs")
                title("u - $(drifter_names[used_ids[idrifter]])")

                subplot(2,1,2)
                plot(time[sel],v[sel]; label = "obs")
                title("v - $(drifter_names[used_ids[idrifter]])")
#            end

            @info "plotting $label"
            subplot(2,1,1)
            plot(time[sel],ui[sel]; label = label)
            subplot(2,1,2)
            plot(time[sel],vi[sel]; label = label)

            mkpath(figdir)
            fig.autofmt_xdate()
            lgd = legend(loc="upper center", bbox_to_anchor=(1.2, 2))
            savefig(joinpath(figdir,"drifter_$(drifter_names[used_ids[idrifter]]).pdf"),bbox_extra_artists=(lgd,), bbox_inches="tight");
            savefig(joinpath(figdir,"drifter_$(drifter_names[used_ids[idrifter]]).png"),bbox_extra_artists=(lgd,), bbox_inches="tight");
        end
    end

    RMSstat = sqrt(nanrms(ui,u)^2 +  nanrms(vi,v)^2)
    @show sum(isfinite.(ui))

    #@show size(ui)
    return (ui,vi),(RMSstat_drifter_u,RMSstat_drifter_v),RMSstat_drifter,RMSstat
end


function validate_drifter(
    fnamerec::AbstractString,
    fnameobs::AbstractString; kwargs...)

    (lon,lat,time),(u,v),ids,drifter_names =
        NCDataset(fnameobs) do ds
            (ds["lon"][:],ds["lat"][:],ds["time"][:]),
            (ds["u"][:],ds["v"][:]),ds["ids"][:],
            ds["drifter_names"][:]
        end


    (ui,vi),(RMSstat_drifter_u,RMSstat_drifter_v),RMSstat_drifter,RMSstat = validate_drifter(
        fnamerec, (lon,lat,time),(u,v),ids;
        drifter_names = drifter_names,
        kwargs...)

    fn = joinpath(dirname(fnamerec),"drifter_validation.nc")

    @info "RMSstat = $RMSstat ($fnamerec)"
    NCDataset(fn,"c") do ds
        defVar(ds,"ui",ui,("obs",))
        defVar(ds,"vi",ui,("obs",))
        defVar(ds,"RMSstat_drifter_u",RMSstat_drifter_u,("drifter",))
        defVar(ds,"RMSstat_drifter_v",RMSstat_drifter_v,("drifter",))
        defVar(ds,"RMSstat_drifter",RMSstat_drifter,("drifter",))
        ds.attrib["RMSstat"] = RMSstat
    end

    return (ui,vi),(RMSstat_drifter_u,RMSstat_drifter_v),RMSstat_drifter,RMSstat
end
