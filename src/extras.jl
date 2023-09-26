using PyPlot
using DIVAnd
using FileWatching
using ForecastVerification
using JLD2
using LaTeXTabulars
using OceanPlot
using Printf
using NCDatasets

include("linreg.jl")

#=
function PyObject(a::Array{Union{T,Missing},N}) where {T,N}
  numpy_ma = pyimport("numpy").ma
  pycall(numpy_ma.array, Any, coalesce.(a,zero(T)), mask=ismissing.(a))
end
=#

function PyObject(a::Array{Union{Missing},N}) where {N}
  numpy_ma = pyimport("numpy").ma
  pycall(numpy_ma.array, Any, coalesce.(a,1.), mask=ismissing.(a))
end

include("uread.jl")
export gread


const fname_all_cv = expanduser("~/Data/Med/AVHRR/Data/avhrr_sub_add_clouds.nc")
const fname_orig = expanduser("~/Data/Med/AVHRR/Data/avhrr_sub3.nc")

const AVHRR_case = (
    fname_orig = expanduser("~/Data/Med/AVHRR/Data/avhrr_sub3.nc"),
    fname_cv = expanduser("~/Data/Med/AVHRR/Data/avhrr_sub_add_clouds.nc"),
    varname = "SST",
    ntest = "last-50",
    dineofrec = expanduser("~/Data/Med/AVHRR/Data/sst_all_cloudsadded.rec"),
    ylim_rms = (0.33,0.8),
)



function ncmean(ds,varname)
    if haskey(ds,"batch_m_rec")
        # old
        return "batch_m_rec"
    elseif haskey(ds,"mean_rec")
        # new
        return "mean_rec"
    else
        # for dineof
        return varname
    end
end

function ncsigma(ds,varname)
    if haskey(ds,"batch_m_sigma")
        # old
        return "batch_m_sigma"
    elseif haskey(ds,"sigma_rec")
        return "sigma_rec"
    elseif haskey(ds,varname * "_error")
        # new
        return varname * "_error"
    else
        return nothing
    end
end

plotres(fname) = plotres(fname,fname_orig)

function plotres(fname::AbstractString,fname_all; dineofrec = [
    expanduser("~/Data/Med/AVHRR/Data/sst_all_cloudsadded.rec"),
    expanduser("~/Data/Med/AVHRR/Data/sst_all_cloudsadded_2008_2009.rec")], offset = -3)
    function myplot(x, t, cl = extrema(skipmissing(x)), cmap="jet")
        pcolor(lon,lat,x'; cmap=cmap)
        set_aspect_ratio()
        title(t)
        clim(cl)
        colorbar()
        plotmap()
        xlim(extrema(lon))
        ylim(extrema(lat))

        return cl
    end

    case = AVHRR_case

    #fname = "data-2019-03-21T101900.nc";

    ntest = 50

    ds = Dataset(fname_all);
    time = ds["time"][end-ntest+1:end]
    close(ds)

    lon,lat,batch_m_true,batch_m_in,batch_m_rec,batch_sigma_rec,mask = loadbatch(case,fname)

    dineof_batch_m_rec = Vector{Any}(undef,length(dineofrec))
    for l = 1:length(dineofrec)
        dineof_batch_m_rec[l] = DINCAE_utils.gread(dineofrec[l],missing)[:,end:-1:1,end-ntest+1:end];
    end

    i = 1
    @show time[i]
    if offset != 0
        @info "offset: $offset"
    end
    batch_m_true .+= offset
    batch_m_in .+= offset
    batch_m_rec .+= offset
    dineof_batch_m_rec[1] .+= offset
    dineof_batch_m_rec[2] .+= offset

    fig = figure(figsize=(10,7))

    for i = 1:length(time)
#    for i = 1:1
        clf()

        cl = extrema(vcat(collect(skipmissing(batch_m_true[:,:,i])),
                          collect(skipmissing(dineof_batch_m_rec[1][:,:,i]))
                          ))

        fig.suptitle("Date: $(Dates.format(time[i],"yyyy-mm-dd"))", fontsize=16)
        fig.text(0.8, 0.05, "GHER, ULiege", fontsize=12)

        #subplot(2,3,1); myplot(batch_m_true[:,:,i],"original $(case.varname) ",cl)
        subplot(2,3,1); myplot(batch_m_true[:,:,i],"(a) original $(case.varname)",cl)
        subplot(2,3,2); myplot(batch_m_in[:,:,i],"(b) original $(case.varname) with added clouds",cl)
        subplot(2,3,3); myplot(batch_m_rec[:,:,i],"(c) DINCAE $(case.varname)",cl)

        diff = batch_m_rec[:,:,i] - batch_m_true[:,:,i]
        mdiff = maximum(abs.(skipmissing(diff)))

        cl_diff = (0, maximum(abs.(skipmissing(batch_sigma_rec[:,:,i]))))
        cl_diff = nothing
        #cl = (-mdiff,mdiff)
        #subplot(2,3,4); myplot(batch_sigma_rec[:,:,i],"exp. error of rec. (std. dev.)",cl_diff,"Reds")
        subplot(2,3,4); myplot(batch_sigma_rec[:,:,i],"(d) exp. error of rec. (std. dev.)",cl_diff,"jet")
        #subplot(2,3,4); myplot(batch_sigma_rec[:,:,i],"exp. error of rec. (std. dev.)",cl_diff,"bwr")

        #subplot(2,3,5); myplot(diff,"true - reconstruction",(-mdiff,mdiff),"bwr")


        subplot(2,3,5); myplot(dineof_batch_m_rec[1][:,:,i],"(e) DINEOF (using all data)",cl)
        subplot(2,3,6); myplot(dineof_batch_m_rec[2][:,:,i],"(f) DINEOF (using 2008-2009)",cl)

        #savefig(replace(fname,".nc" => @sprintf("-%03d-new.png",i)))
        savefig(replace(fname,".nc" => @sprintf("-%03d-new2-300.png",i)),dpi=300)
    end
end


loadall(fname) = loadall(AVHRR_case,fname)

function loadall(case,fname)
    ds = Dataset(fname);

    mean_varname = ncmean(ds,case.varname)
    batch_m_rec = ds[mean_varname][:,:,:]

    # add mean
    if haskey(ds,"batch_m_rec")
        # old files
        meandata = ds["meandata"][:][:,:,1:1]
        batch_m_rec = batch_m_rec .+ meandata
    end

    close(ds)

    return batch_m_rec
end


function fixdineof(case,fname)
    lon,lat,time,time_units =
        Dataset(case.fname_orig) do ds
            ds["lon"][:],ds["lat"][:],ds["time"][:],ds["time"].attrib["units"]
        end


    Dataset(fname,"a") do ds
#=    renameDim(ds,"dim001","lon");
    renameDim(ds,"dim002","lat");
    renameDim(ds,"dim004","time");
    renameVar(ds,"sst_t","mean_rec")
=#
        defVar(ds,"lon",lon,("lon",))
        defVar(ds,"lat",lat,("lat",))
        defVar(ds,"time",time,("time",), attrib = ["units" => time_units])
    end
    return nothing
end

recavg(fnames,fnameavg) = recavg(AVHRR_case,fnames,fnameavg)

function recavg(case,fnames,fnameavg)
    cp(fnames[1],fnameavg, force=true)

    ds = Dataset(fnames[1])

    mean_varname = ncmean(ds,case.varname)
    sigma_varname = ncsigma(ds,case.varname)
    sz = size(ds[mean_varname])

    close(ds)

    batch_m_rec = zeros(sz)
    batch_sigma_rec = zeros(sz)

    for i = 1:length(fnames)
        @show i
        ds = Dataset(fnames[i])
        batch_m_rec += ds[mean_varname][:]
        batch_sigma_rec += ds[sigma_varname][:]
        close(ds)
    end

    batch_sigma_rec /= length(fnames)
    batch_m_rec /= length(fnames)

    Dataset(fnameavg,"a") do ds
        ds[mean_varname][:] = batch_m_rec
        ds[sigma_varname][:] = batch_sigma_rec
    end
end


function recmedian(fnames,fnamemedian)
    cp(fnames[1],fnamemedian, force=true);

    ds = Dataset.(fnames);

    Dataset(fnamemedian,"a") do dsm
        chunk = 100
        tmp = zeros(Union{Missing,Float32},dsm.dim["lon"],dsm.dim["lat"],length(fnames),chunk);
        ntime = dsm.dim["time"]

        for varname = ["batch_m_rec","batch_sigma_rec"]

            for n = 1:chunk:ntime
                n2 = min(n+chunk-1,ntime)
                if n % 1 == 0
                    @info "n = $n"
                end

                for i = 1:length(fnames)
                    tmp[:,:,i,1:(n2-n+1)] = ds[i][varname][:,:,n:n2];
                end

                dsm[varname][:,:,n:n2] = replace(median(replace(tmp,missing => NaN),dims = 3),NaN => missing);
            end
        end
    end
    close.(ds)
end

function loadbatchavg(case,fnames)
    lon,lat,batch_m_true,batch_m_in,batch_m_rec,batch_sigma_rec,mask = DINCAE_utils.loadbatch(case,fnames[1])

    for i = 2:length(fnames)
        lon,lat,batch_m_true_,batch_m_in_,batch_m_rec_,batch_sigma_rec_,mask = DINCAE_utils.loadbatch(case,fnames[i])

        batch_m_rec += batch_m_rec_
        batch_sigma_rec += batch_sigma_rec_
    end
    batch_sigma_rec /= length(fnames)
    batch_m_rec /= length(fnames)

    return lon,lat,batch_m_true,batch_m_in,batch_m_rec,batch_sigma_rec,mask
end

errstat(fname::AbstractString) = errstat(stdout,fname)



function summarydineof(dineofrec)
    case = AVHRR_case
    fnamesummary = dineofrec * ".json"

    if !isfile(fnamesummary)
        ntest = 50

        batch_m_rec = DINCAE_utils.gread(dineofrec,missing)[:,end:-1:1,end-ntest+1:end];
        ds = Dataset(case.fname_orig)
        batch_m_true = ds[case.varname][:,:,end-ntest+1:end];
        mask = ds["mask"][:];
        close(ds)

        ds = Dataset(case.fname_cv)
        batch_m_in = ds[case.varname][:,:,end-ntest+1:end];
        close(ds)

        mm = ismissing.(batch_m_in) .& .!ismissing.(batch_m_true) .& reshape(mask .== 1,(size(mask,1),size(mask,2),1));

        m_true = batch_m_true[mm]
        m_rec = batch_m_rec[mm]
        RMS = rms(m_true,m_rec)

        summary = Dict(
            "cvrms" => RMS,
            "std_true" => std(m_true),
            "std_rec" => std(m_rec),
            "cor" => cor(m_true,m_rec),
            "cvcrms" => crms(m_true,m_rec),
            "number" => sum(mm)
        )


        open(fnamesummary,"w") do f
            JSON.print(f, summary)
        end
        return summary
    else
        return JSON.parse(read(fnamesummary,String))
    end
end

function rmsdineof(dineofrec)
    case = AVHRR_case
    ntest = 50

    batch_m_rec = DINCAE_utils.gread(dineofrec,missing)[:,end:-1:1,end-ntest+1:end];
    ds = Dataset(case.fname_orig)
    batch_m_true = ds[case.varname][:,:,end-ntest+1:end];
    mask = ds["mask"][:];
    close(ds)

    ds = Dataset(case.fname_cv)
    batch_m_in = ds[case.varname][:,:,end-ntest+1:end];
    close(ds)

    mm = ismissing.(batch_m_in) .& .!ismissing.(batch_m_true) .& reshape(mask .== 1,(size(mask,1),size(mask,2),1));

    m_true = batch_m_true[mm]
    m_rec = batch_m_rec[mm]
    RMS = rms(m_true,m_rec)

    return RMS
end

function RMS_iteration(case; fnameavg = nothing, yl = case.ylim_rms)
    clf();
    RMScv = [DINCAE_utils.cvrms(case,fn) for fn in sort(glob("data*T*.nc"))]
    plot(10*(1:length(RMScv)),RMScv,label="DINCAE");
    xr = xlim()

    if haskey(case,:dineofrec)
        RMSdineof = DINCAE_utils.rmsdineof(case.dineofrec)
        hlines(RMSdineof,xr[1],xr[2],linestyle="--",colors="c",label="DINEOF");
    end

    if fnameavg != nothing
        RMScvavg = DINCAE_utils.cvrms(case,fnameavg)
        hlines(RMScvavg,xr[1],xr[2],linestyle="--",colors="r",label="DINCAE avg. rec.");
    end

    legend()
    xlabel("Iteration (epoch)");
    ylabel("RMS (cross-validation)");
    if yl != nothing
        ylim(yl)
    end
    savefig("RMS_iterations.png");
    savefig("RMS_iterations.pdf");
end

function compareall(; dineofrec = expanduser("~/Data/Med/AVHRR/Data/sst_all_cloudsadded.rec"))
    fnameavg = "data-avg.nc"
    fnamemedian = "data-median.nc"
    dineofrec_ = expanduser("~/Data/Med/AVHRR/Data/sst_all_cloudsadded_2008_2009.rec")
    #fnameavg = "data-median.nc"
    DINCAE_utils.compare([
        "DINEOF",
        "DINEOF (2008-2009)",
        "DINCAE (no skip connections)",
        "DINCAE (2 skip connections)",
        "DINCAE (all skip connections)",
        "DINCAE (all skip connections - median)",
        "DINCAE (all skip connections and wider layers)",
        "DINCAE (all skip connections and narrower layers)",
        "DINCAE (all skip connections and less layers)",
        "DINCAE (all skip connections and more layers)",
        "DINCAE (all skip connections and average pooling)",
    ],
                         [
                             dineofrec,
                             dineofrec_,
                             joinpath("Fig-jitter-noskip",fnameavg),
                             joinpath("Fig-jitter",fnameavg),
                             joinpath("Fig-jitter-more-skip",fnameavg),
                             joinpath("Fig-jitter-more-skip",fnamemedian),
                             joinpath("Fig-jitter-more-skip-wider-layers",fnameavg),
                             joinpath("Fig-jitter-more-skip-narrow-layers",fnameavg),
                             joinpath("Fig-jitter-more-skip-less-layers",fnameavg),
                             joinpath("Fig-jitter-more-skip-more-layers",fnameavg),
                             joinpath("Fig-jitter-more-skip-avg-pool",fnameavg),
                             #joinpath("Fig-jitter-more-skip-more-layers","Old",fnameavg),
                          ])
end

function compare(cases,fnames)
    case = AVHRR_case

    fnameavg = fnames[end]
    SSTc = DINCAE_utils.loadall(fnameavg);
    dineofrec = expanduser("~/Data/Med/AVHRR/Data/sst_all_cloudsadded.rec")
    SSTd = reverse(DINCAE_utils.gread(dineofrec,missing),dims = 2);
    SSTorig = Dataset(case.fname_orig)[case.varname][:];
    SSTadd = Dataset(case.fname_cv)[case.varname][:];


    function comparecase(name,fnameavg)
        if endswith(fnameavg,".nc")
            SST = DINCAE_utils.loadall(fnameavg);
        else
            SST = reverse(DINCAE_utils.gread(fnameavg,missing),dims = 2);
        end

        n = (size(SSTorig,3) - size(SST,3)+1) :size(SSTorig,3)
        m = ismissing.(SSTadd[:,:,n]) .& (.!ismissing.(SSTorig[:,:,n])) .& .!ismissing.(SSTd[:,:,n]) .& .!ismissing.(SST);

        @info "RMS $name"
        @info "CV ForecastVerification"
        cv_dincae = ForecastVerification.summary(SST[m],SSTorig[:,:,n][m])

        m = .!ismissing.(SSTadd[:,:,n]) .& (.!ismissing.(SSTorig[:,:,n])) .& .!ismissing.(SSTd[:,:,n]) .& .!ismissing.(SST);
        @info "non-CV ForecastVerification"
        noncv_dincae = ForecastVerification.summary(SST[m],SSTorig[:,:,n][m])
        return [name,cv_dincae,noncv_dincae]
    end

    @show cases,fnames
    fmt(x) = @sprintf("%2.4f",x)
    #data = [["DINCAE",cv_dincae,noncv_dincae],
    #        ["DINEOF",cv_dineof,noncv_dineof]]
    #data = comparecase.(["DINCAE","DINEOF"],[dineofrec,fnameavg])
    data = comparecase.(cases,fnames) # ["DINCAE","DINEOF"],[dineofrec,fnameavg])

    tabname = "table-validation2.tex"
    latex_tabular(tabname,
              Tabular("lrrrrrr"),
                  [Rule(:top),
                   ["",MultiColumn(3, :c, "CV data"),MultiColumn(3, :c, "non-CV data")],
                   ["","RMS", "CRMS", "bias","RMS", "CRMS", "bias"],
                   Rule(:mid),
                   map(r -> [lpad(r[1],20),fmt(r[2].rms),fmt(r[2].crms),fmt(r[2].bias),fmt(r[3].rms),fmt(r[3].crms),fmt(r[3].bias)],data)...,
                   Rule(:bottom)])

    println.(read(tabname,String))
end


function spectrum(fnameavg,dineofrec)
    SST = nomissing(loadall(fnameavg),0)

    SSTd = gread(dineofrec,0.)

    krad,myradspec = ForecastVerification.rad_power_spectrum(4_000,4_000,SST)
    kradd,myradspecd = ForecastVerification.rad_power_spectrum(4_000,4_000,SSTd)

    loglog(krad,myradspec,"b-",label="DINCAE")
    loglog(kradd,myradspecd,"r-",label="DINEOF")
    legend();
    savefig("spectrum.png");
end

function taylorplot(fname; dineofrec = expanduser("~/Data/Med/AVHRR/Data/sst_all_cloudsadded.rec"))

    d = DINCAE_utils.summarydineof(dineofrec)
    dc = summary(fname)
    ForecastVerification.taylorplot(
        [d["cvcrms"],dc["cvcrms"]],
        [d["cor"],dc["cor"]],
        [d["std_rec"],dc["std_rec"]],d["std_true"]; labels=["DINEOF","DINCAE"])

    savefig("taylorplot.png")
end


function plotregion()
    figure(figsize=(6,5))
    ds = Dataset(DINCAE_utils.fname_orig);
    lon = ds["lon"][:];
    lat = ds["lat"][:];
    close(ds)

    bathname = expanduser("~/Data/DivaData/Global/gebco_30sec_4.nc");
    bx,by,b = DIVAnd.extract_bath(bathname,true,-6:20,30:45);
    b[b .< 0] .= NaN;
    pcolor(bx,by,b'; cmap="jet");
    set_aspect_ratio();
    colorbar(orientation="horizontal");
    #colorbar()
    lonr = extrema(lon)
    latr = extrema(lat)
    plot(lonr[[1,2,2,1,1]],latr[[1,1,2,2,1]],"r-")
    #xlim(-6,20)
    ylim(30,45.5)
    plotmap()
    savefig("plotregion.png",dpi=300);
end

function seasonalaverage(SST,time; DT = 30)
    cycle_len = 365
    half_DT = DT/2
    half_cycle_len = cycle_len/2
    doy = Dates.dayofyear.(time);
    mSST = zeros(eltype(SST),size(SST,1),size(SST,2),maximum(doy))

    Threads.@threads for j = 1:size(SST,2)
        for i = 1:size(SST,1)
            for nn = 1:size(mSST,3)
            #for nn = 1:10
                count = 0

                for n = 1:length(time)
                    #if isfinite(SST[i,j,n])
                    if !ismissing(SST[i,j,n])
                        if abs( mod( doy[n] - nn + half_cycle_len, cycle_len) - half_cycle_len) <= half_DT
                            mSST[i,j,nn] += SST[i,j,n]
                            count += 1
                        end
                    end
                end
                mSST[i,j,nn] /= count

#                #for i = 1:2
#                sel .= abs.(  mod.( doy  .- doy[nn] .+ 365/2, 365) .- 365/2) .<= DT÷2;
#                mSST[i,j,nn] = mean(@view SST[i,j,sel])
            end
        end
    end
    return mSST

end


function remove_seasonal_cycle(SST,SSTtime; DT = 30)
    doy = Dates.dayofyear.(SSTtime);
    mSST2 = DINCAE_utils.seasonalaverage(SST,SSTtime; DT = DT);

    SSTa = similar(SST);
    for n = 1:size(SST,3)
        SSTa[:,:,n] = SST[:,:,n] - mSST2[:,:,doy[n]]
    end
    return SSTa
end

function std_around_seasonalaverage(SST,SSTtime)
    doy = Dates.dayofyear.(SSTtime);
    mSST2 = DINCAE_utils.seasonalaverage(SST,SSTtime);

    SSTa = similar(SST);
    for n = 1:size(SST,3)
        SSTa[:,:,n] = SST[:,:,n] - mSST2[:,:,doy[n]]
    end

    count = sum(.!ismissing.(SSTa), dims = 3)
    SSTa[ismissing.(SST)] .= 0

    SST_std = sqrt.(sum(SSTa.^2,dims = 3) ./ count)[:,:,1];
    return SST_std
end


function plot_std_around_seasonalaverage(lon,lat,SST_std; vmin = 0.7, vmax = 1.4)
    ima = pcolor(lon,lat,SST_std',cmap= "jet",vmin = vmin, vmax = vmax);
    set_aspect_ratio();
    #colorbar(orientation="horizontal")
    plotmap()
    xlim(lon[1],lon[end])
    ylim(lat[end],lat[1])
    return ima
end


function compare_std_around_seasonalaverage(;
                                            dineofrec = expanduser("~/Data/Med/AVHRR/Data/sst_all_cloudsadded.rec"),
                                            fnameavg = "data-avg.nc",
                                            fname_orig = "/home/abarth/Data/Med/AVHRR/Data/avhrr_sub3.nc"
                                            )

    case = AVHRR_case
    SSTdineof = reverse(DINCAE_utils.gread(dineofrec,missing),dims = 2)
    SSTdincae = DINCAE_utils.loadall(fnameavg);

    ds = Dataset(fname_orig)
    SSTorig = ds[case.varname][:,:,:];
    lon = nomissing(ds["lon"][:]);
    lat = nomissing(ds["lat"][:]);
    SSTtime = nomissing(ds["time"][:]);
    close(ds)

    figure(figsize=(10,6))
    SSTdineof_std = @time DINCAE_utils.std_around_seasonalaverage(SSTdineof,SSTtime);
    SSTdincae_std = @time DINCAE_utils.std_around_seasonalaverage(SSTdincae,SSTtime);
    SSTorig_std = @time DINCAE_utils.std_around_seasonalaverage(SSTorig,SSTtime);

    subplot(1,3,1)
    DINCAE_utils.plot_std_around_seasonalaverage(lon,lat,SSTorig_std)
    title("original $(case.varname) std. dev.")

    subplot(1,3,2)
    DINCAE_utils.plot_std_around_seasonalaverage(lon,lat,SSTdincae_std)
    title("DINCAE $(case.varname) std. dev.")

    subplot(1,3,3)
    ima = DINCAE_utils.plot_std_around_seasonalaverage(lon,lat,SSTdineof_std)
    title("DINEOF $(case.varname) std. dev.")

    subplots_adjust(bottom = 0.1);
    cbar_ax = gcf().add_axes([0.3, 0.15, 0.4, 0.025])
    colorbar(ima, cax=cbar_ax, orientation = "horizontal");

    savefig("compare_std_around_seasonalaverage.png",dpi=300);
end


function post_process()
    fnames = sort(glob("*T*.nc"))[20:100]
    fnameavg = "data-avg.nc"
    #DINCAE_utils.recavg(fnames,fnameavg)
    DINCAE_utils.errstat(fnameavg)
    DINCAE_utils.plotres(fnameavg)
    compare_std_around_seasonalaverage();
end

post_process_sst() = post_process(cases("sst_t"))
post_process_chlor_a() = post_process(cases("chlor_a"))
post_process(case::String) = post_process(cases(case))
post_process2(case::String) = post_process2(cases(case))

function post_process(case)
    fnames = sort(glob("*T*.nc"))[20:end]
    fnameavg = "data-avg.nc"
    recavg(fnames,fnameavg)
    RMS_iteration(case; fnameavg = fnameavg, yl = case.ylim_rms)
    errstat(stdout,case,fnames[end])
end

function post_process2(case)
    fnames = sort(glob("*T*.nc"))[10:100]
    fnameavg = "data-avg2.nc"
    recavg(fnames,fnameavg)
    RMS_iteration(case; fnameavg = fnameavg, yl = case.ylim_rms)
end

watch_post_process(case::String; path = pwd()) = watch_post_process(cases(case); path = pwd())

function watch_post_process(case; path = path)
    @info "watching $path"
    while true
        fname,status = FileWatching.watch_folder(path)
        # make sure file is saved
        sleep(20)
        @show DINCAE_utils.cvrms.(Ref(case),sort(glob("*T*.nc")));
        @show Dates.now()
    end
end


function cases(p)
    basedir = expanduser("~/Data/DINCAE-multivariate/Adriatic2/")

            #fname_rec: sort(glob("*T*.nc"))[end]
    if p == "chlor_a"
        return (
            fname_orig = joinpath(basedir,"color_revlat.nc"),
            fname_cv = joinpath(basedir,"color_revlat_add_clouds.nc"),
            varname = "chlor_a",
            transfun = (log,exp),
            clim = (0.1,50),
            ylim_rms = (0.2, 2),
            ntest = "all",
        )
    elseif p == "chlor_a_log"
        return (
            fname_orig = joinpath(basedir,"color_revlat_log.nc"),
            fname_cv = joinpath(basedir,"color_revlat_log_add_clouds.nc"),
            varname = "chlor_a",
            transfun = (identity,identity),
            clim = (-2.,4),
            ylim_rms = nothing,
            ntest = "all"
        )
    elseif p == "sst_t"
        return (
            fname_orig = joinpath(basedir,"modis_sst_revlat.nc"),
            fname_cv = joinpath(basedir,"modis_sst_revlat_add_clouds.nc"),
            varname = "sst_t",
            transfun = (identity,identity),
            clim = (10,28),
            ylim_rms = (0.6,2.4),
            ntest = "all"
        )
    elseif p == "uwnd"
        return (
            fname_orig = joinpath(basedir,"CCMP_Wind_Analysis_Adriatic_revlat.nc"),
#            fname_cv = joinpath(basedir,"CCMP_Wind_Analysis_Adriatic_revlat.nc"),
            varname = "uwnd",
            transfun = (identity,identity),
            clim = (10,28),
            ylim_rms = (0.75,2.4),
            ntest = "all"
        )
    elseif p == "vwnd"
        return (
            fname_orig = joinpath(basedir,"CCMP_Wind_Analysis_Adriatic_revlat.nc"),
#            fname_cv = joinpath(basedir,"CCMP_Wind_Analysis_Adriatic_revlat.nc"),
            varname = "vwnd",
            transfun = (identity,identity),
            clim = (10,28),
            ylim_rms = (0.75,2.4),
            ntest = "all"
        )
    end
end
function prep_mv()
    Random.seed!(1234)
    #varname = "chlor_a"
    varname = "sst_t"

    case = cases(varname)

    if varname == "sst_t"
        #prep_sst()
    end

    fname_cv = addcvpoint(case.fname_orig,case.varname)


end


function loadvalid(fname,varname,qualname,varrange,qualrange)
    Dataset(fname,"r") do ds
        var = ds[varname][:,:,:]
        qual = ds[qualname][:,:,:]

        good = (qualrange[1] .<= qual .<= qualrange[2]) .& .!ismissing.(qual) .& (
                varrange[1] .<= var .<= varrange[2])

        var[.!good] .= missing
        return var[:,:,1]
    end
end



function prep_sst()
    fname = "modis_sst_revlat.nc"

    Dataset(fname,"a") do ds
        sst = ds["sst"][:,:,:]
        qual = ds["qual"][:,:,:]
        sst_t = copy(sst)
        sst_t[(qual .> 3) .& .!ismissing.(qual)] .= missing
        sst_t[.!ismissing.(sst) .& (sst_t .> 40)] .= missing

        if haskey(ds,"sst_t")
            ds["sst_t"][:] = sst_t
        else
            defVar(ds,"sst_t", sst_t, ("lon", "lat", "time"), fillvalue = Float32(-9999.0))
        end
    end
    return nothing
end

function prep_chlor_a_log()
    case = DINCAE_utils.cases("chlor_a")
    case2 = DINCAE_utils.cases("chlor_a_log")
    cp(case.fname_orig,case2.fname_orig)
    cp(case.fname_cv,case2.fname_cv)

    ds = Dataset(case2.fname_orig,"a");
    ds[case2.varname][:,:,:] = log.(ds[case2.varname][:,:,:]);
    close(ds)

    ds = Dataset(case2.fname_cv,"a");
    ds[case2.varname][:,:,:] = log.(ds[case2.varname][:,:,:]);
    close(ds)
end

function addmask!(data,mask)
    for n = 1:size(data,3)
        tmp = data[:,:,n]
        tmp[.!mask] .= missing
        data[:,:,n] = tmp
    end
end

function prep_mv2()
    mincvfrac = 0.5
    split = [("train",0.7),("dev",0.2),("test",0.1)]
    #split = [("test",0.1)]

    data = [
       (filename = "/home/abarth/Data/DINCAE-multivariate/Adriatic2/modis_sst_revlat.nc",
        varname = "sst_t"),
       (filename = "/home/abarth/Data/DINCAE-multivariate/Adriatic2/color_revlat_log.nc",
        varname = "chlor_a"),]

    for d in data
        add_mask(d.filename,d.varname)

        newfnames = DINCAE_utils.splitdata(d.filename,split)

        addcvpoint(newfnames[2],d.varname, mincvfrac = mincvfrac)
        addcvpoint(newfnames[3],d.varname, mincvfrac = mincvfrac)
    end
end



function corp(sst,wnd)
    sel = findall(.!ismissing.(sst))
    a,b,R2 = linreg(wnd[sel,:],sst[sel])
    #return cor(sst[sel],wnd[sel])
    return R2
end



# jmax iterations of a Laplacian smoother

#
# for x = 0, expept x[500] = 1
# the smoothed solution is
#
# 1/(sqrt(pi * jmax)) * exp.(  -(x2 .- 500).^2 / ( jmax ))
#
# the standard deviation of this gaussian is
# sqrt(jmax/2)
#
# for a given standard deviation σ, jmax is 2σ²


function smooth!(x,jmax)
    imax = length(x)

    for j = 1:jmax
        xm = x[1]

        for i = 2:imax-1
            xf = 0.5 * x[i] + 0.25 * (x[i+1] + xm)
            xm = x[i]
            x[i] = xf
        end
    end

end

function lagcorr(sst,wnd,lag)
    n = size(sst,3)
    n1 = max(1,1+lag) : min(n+lag,n);
    n2 = max(1,1-lag) : min(n-lag,n);
    return corp((@view sst[:,:,n1]),(@view wnd[:,:,n2,:]))
end

function testlags(sst,wnd,lags)

    CC = zeros(length(lags))
    Threads.@threads for l = 1:length(lags)
        CC[l] = lagcorr(sst,wnd,lags[l]);
    end
    return CC
end

function testlagssmooth(sst,wnd_,lags,smoothing_jmax)
    wnd = copy(wnd_)
    jmax = 0
    CC = zeros(length(lags),length(smoothing_jmax))

    for k = 1:length(smoothing_jmax)
        @info("smooth for $(smoothing_jmax[k]-jmax)")

        for n = 1:size(wnd,4)
            Threads.@threads for j = 1:size(wnd,2)
                for i = 1:size(wnd,1)
                    smooth!((@view wnd[i,j,:,n]),smoothing_jmax[k]-jmax)
                end
            end
        end
        jmax = smoothing_jmax[k]

        #CC = testlags(sst,wnd,kags)
        CC[:,k] = testlags(sst,wnd,lags)
        @info("maximum lag $(maximum(abs.(CC[:,k])))")
    end

    return CC
end


function testlags()
    res = testlags.([:speed,1,2,3])
    return res
end

function testlags(windtype)
    sst_case = DINCAE_utils.cases("sst_t")
    uwnd_case = DINCAE_utils.cases("uwnd")
    vwnd_case = DINCAE_utils.cases("vwnd")

    uwnd =
        Dataset(uwnd_case.fname_orig) do ds
            ds[uwnd_case.varname][:,:,:];
        end

    vwnd =
        Dataset(vwnd_case.fname_orig) do ds
            ds[vwnd_case.varname][:,:,:];
        end

    sst,datatime =
        Dataset(sst_case.fname_orig) do ds
            ds[sst_case.varname][:,:,:],ds["time"][:]
        end

    #wnd = uwnd + 1im * vwnd
    #wnd = abs.(uwnd + 1im * vwnd)
    #wnd = fun.(uwnd,vwnd)
    DT = 3
    ssta = DINCAE_utils.remove_seasonal_cycle(sst,datatime; DT = DT)
    #uwnda = DINCAE_utils.remove_seasonal_cycle(uwnd,datatime; DT = DT)
    #vwnda = DINCAE_utils.remove_seasonal_cycle(vwnd,datatime; DT = DT)

    #wnda = cat(uwnda,vwnda,dims=Val(4))
    speed = sqrt.(uwnd.^2 + vwnd.^2)

    wnd =
        if windtype == :speed
            speed
        else
            nexp = windtype::Int
            cat((speed.^nexp) .* uwnd, (speed.^nexp) .* vwnd,dims=Val(4))
        end

    lags = -10:30;
    lags = -10:30;
    lags = 30:60;
    lags = -10:60;
    #lags = -1:1;

    σ_time = 0:5:50

    lags = -5:10;
    σ_time = 0:1:10

    lags = -2:10;
    σ_time = 0:0.5:5

    lags = -2:15;
    σ_time = 0:0.5:15

    smoothing_jmax = round.(Int,2 * σ_time.^2)
    @show smoothing_jmax

    smoothing_jmax = 0:2:50
    σ_time = sqrt.(smoothing_jmax/2)
    @show σ_time

    CC = DINCAE_utils.testlagssmooth(ssta,wnd,lags,smoothing_jmax)

    #plot(lags,CC); ylabel("correlation coefficient"); xlabel("time shift in days")
    #pcolor(lags,smoothing_jmax,CC'); ylabel("correlation coefficient"); xlabel("time shift in days")
    #ylabel("smoothing")
    return lags,smoothing_jmax,σ_time,CC

    #return lags,CC
end


function prep_wind_speed()
    #JLD2.@load "DINCAE_utils.testlags_speed123.jld2"
    res =
        JLD2.jldopen("DINCAE_utils.testlags_speed123.jld2") do file
            file["res"]
        end

    for (lags,smoothing_jmax,σ_time,CC) in res
        i,j = Tuple(findmax(CC)[2])
        @show findmax(CC)[1],size(CC),length(σ_time),lags[i],σ_time[j]
    end

    l = 1
    lags,smoothing_jmax,σ_time,CC = res[l]
    i,j = Tuple(findmax(CC)[2])

    @show findmax(CC)[1],size(CC),length(σ_time),lags[i],σ_time[j]
    uwnd_case = DINCAE_utils.cases("uwnd")
    vwnd_case = DINCAE_utils.cases("vwnd")

    uwnd,vwnd,lon,lat,datatime =
        Dataset(uwnd_case.fname_orig) do ds
            (ds[uwnd_case.varname][:,:,:],
             ds[vwnd_case.varname][:,:,:],
             ds["lon"][:],
             ds["lat"][:],
             ds["time"][:])
        end

    speed = sqrt.(uwnd.^2 + vwnd.^2)
    windtype = [:speed,1,2,3][l]

    wnd =
        if windtype == :speed
            speed
        else
            nexp = windtype::Int
            cat((speed.^nexp) .* uwnd, (speed.^nexp) .* vwnd,dims=Val(4))
        end

    lag = lags[i]

    for n_ = 1:size(wnd,4)
        Threads.@threads for j_ = 1:size(wnd,2)
            for i_ = 1:size(wnd,1)
                smooth!((@view wnd[i_,j_,:,n_]),smoothing_jmax[j])
            end
        end
    end

    n = size(wnd,3)
    n1 = max(1,1+lag) : min(n+lag,n);
    n2 = max(1,1-lag) : min(n-lag,n);
    wnd_lag = copy(wnd)
    wnd_lag[:,:,n1,:] = wnd[:,:,n2,:]

    ds = Dataset("CCMP_Wind_Analysis_Adriatic_revlat_speed_filtered.nc","c")
    # Dimensions

    ds.dim["lon"] = size(wnd,1)
    ds.dim["lat"] = size(wnd,2)
    ds.dim["time"] = Inf # unlimited dimension

    # Declare variables

    nclon = defVar(ds,"lon", Float64, ("lon",))
    nclat = defVar(ds,"lat", Float64, ("lat",))
    nctime = defVar(ds,"time", Float64, ("time",))
    nctime.attrib["units"] = "days since 1900-01-01 00:00:00"
    ncwnd = defVar(ds,"wind_speed", Float32, ("lon", "lat", "time"))
    ncwnd.attrib["_FillValue"] = Float32(-9999.0)
    ncwnd.attrib["comments"] = "lag $(lags[i]) filter $(σ_time[j])"


    # Define variables

    nclon[:] = lon
    nclat[:] = lat
    nctime[:] = datatime
    ncwnd[:] = wnd

    close(ds)

end


function compare_std_time_serie(case,fname)
    data = Dataset(fname) do ds
        if haskey(ds,case.varname)
            data = ds[case.varname][:,:,:];
            data[data .== 9999] .= missing;
            data
        else
            ds["mean_rec"][:,:,:];
        end
    end


    variance = [var(skipmissing(data[:,:,n])) for n = 1:size(data,3)]

    return sqrt(mean(variance)),variance
end



"""

Example
```julia

filename = "/home/abarth/tmp/DINCAE.jl/MultiSync-chl-refine-lwr-1/MultiSync-CHL_SNS_part_ini_add_clouds.nc"
varname = "CHL"
DINCAE_utils.stat(filename,varname)

```
"""
function stat(filename,varname)
    NCDataset(filename) do ds
        data = ds[varname][:,:,:]
        mask = ds["mask"][:,:]
        n = length(mask)

        sum_mask = sum(mask)
        println("Sea points:  $(sum_mask) out of $n or $(100*sum_mask/n) %")
        println("Land points: $(n-sum_mask) out of $n or $(100*(n-sum_mask)/n) %")

        nvalid = mapslices(slice -> sum(.!ismissing.(slice) .& (mask .== 1)),
                           data,dims=(1,2));

        nvalid_total = sum(nvalid)
        nsea_total = sum(mask)*size(data,3)
        nvalid_percentage = 100 * nvalid_total / nsea_total

        println("Sea points with measurements $(nvalid_total) out of $(nsea_total) or $nvalid_percentage %")
    end
end



function plotstat(
    filename;
    orientation="vertical", resolution='f',
    grid=1.25, title = "Percentage of valid data")

    ds = NCDataset(filename)
    lon = ds["lon"][:]
    lat = ds["lat"][:]

    count_nomissing = ds["count_nomissing"][:,:]
    ntime = ds.dim["time"]::Int

    clf();
    pcolormesh(lon,lat,100 * count_nomissing'/ntime);
    colorbar(orientation=orientation)
    PyPlot.title(title)

    mlon,mlat,mdata = GeoDatasets.landseamask(resolution=resolution, grid=grid)

    contourf(mlon,mlat,mdata',levels  = [0.5,2],colors = [[.5,.5,.5]])
    xlim(extrema(lon))
    ylim(extrema(lat))
    set_aspect_ratio()
    close(ds)
end



function errstat_time(case,fname_rec)
    fname_orig = case.fname_orig
    fname_cv = case.fname_cv
    varname = case.varname

    ds_rec = NCDataset(fname_rec,"r")
    ds_cv = Dataset(fname_cv);
    ds_orig = Dataset(fname_orig);

    time = ds_orig["time"][:]
    lon = ds_orig["lon"][:]
    lat = ds_orig["lat"][:]
    mask = ds_cv["mask"][:,:] .== 1


    rms_cv = fill(NaN,length(time))
    rms_noncv = fill(NaN,length(time))
    bias_cv = fill(NaN,length(time))
    q10 = fill(NaN,length(time))
    q90 = fill(NaN,length(time))
    min_err = fill(NaN,length(time))
    max_err = fill(NaN,length(time))
    min_abs_err = fill(NaN,length(time))
    max_abs_err = fill(NaN,length(time))
    max_abs_err_lon = fill(NaN,length(time))
    max_abs_err_lat = fill(NaN,length(time))

    rr(x) = replace(x, NaN => missing)
    n = 80
    for n = 1:length(time)
        data_cv = ds_cv[varname][:,:,n]
        data_orig = ds_orig[varname][:,:,n]
        data_rec = ds_rec[varname][:,:,n]
        data_orig[.!mask] .= missing

        difference = data_rec - data_orig
        mask_cv = ismissing.(data_cv) .&  .!ismissing.(data_orig)
        mask_noncv = .!ismissing.(data_cv) .&  .!ismissing.(data_orig)
        sum(mask_cv)
        diff_cv = Float64.(difference[mask_cv])
        diff_noncv = Float64.(difference[mask_noncv])

        if length(diff_cv) > 0
            bias_cv[n] = mean(diff_cv)
            rms_cv[n] = sqrt(mean(diff_cv.^2))
            rms_noncv[n] = sqrt(mean(diff_noncv.^2))

            q10[n],q90[n] = quantile(abs.(diff_cv),(0.1,0.9))
            min_err[n] = minimum(diff_cv)
            max_err[n] = maximum(diff_cv)

            min_abs_err[n] = minimum(abs.(diff_cv))
            max_abs_err[n] = maximum(abs.(diff_cv))

            idx = findmax(coalesce.(abs.(difference),-Inf))[2]
            max_abs_err_lon[n] = lon[idx[1]]
            max_abs_err_lat[n] = lat[idx[2]]
        end
    end

    close(ds_orig)
    close(ds_cv)
    close(ds_rec)

    errstat = (; rms_cv, bias_cv, q10, q90, min_err, max_err, min_abs_err, max_abs_err, max_abs_err_lon, max_abs_err_lat, rms_noncv)

    return (time,errstat)
end



function plotres()
    varname = "chlor_a"
    varname = "sst_t"

    case = cases(varname)
    fname_rec = sort(glob("*T*.nc"))[end]

    plotres(case,fname_rec; transfun = case.transfun, clim = case.clim)
end


_noncloudstat(x,mask) = sum(isfinite.(x) .& mask) / (sum(mask) * size(x,3))
_cloudstat(x,mask) = 1 - _noncloudstat(x,mask)

function cloudstat(case::NamedTuple)

    SST_orig,mask = Dataset(case.fname_orig) do ds
        nomissing(ds[case.varname][:,:,:],NaN), ds["mask"][:,:] .== 1
    end;
    SST_cv = Dataset(case.fname_cv) do ds
        nomissing(ds[case.varname][:,:,:],NaN);
    end;

    ntest = 50

    @show 100 * _cloudstat(SST_orig,mask)
    @show 100 * _cloudstat(SST_cv,mask)

    @show 100 * _cloudstat(SST_orig[:,:,end-ntest+1:end],mask)
    @show 100 * _cloudstat(SST_cv[:,:,end-ntest+1:end],mask)
end
