
function splitindices(totlength,split = [("train",0.7),("dev",0.2),("test",0.1)])
    ranges = []
    i = 0

    for (name,fraction)  in split
        len = round(Int,totlength * fraction)
        range =  (i+1) : min(i+len, totlength)
        i += len
        push!(ranges,(;range, name))
        println("$name: indices: $range; fraction: $(length(range)/totlength)")
    end

    return ranges
end



"""
    DINCAE_utils.splitdata(fname,split = [("train",0.7),("dev",0.2),("test",0.1)])

Split the NetCDF file `fname` into training, developpement and test dataset
with the specified fraction.
"""
function splitdata(fname,split = [("train",0.7),("dev",0.2),("test",0.1)])
    ds = NCDataset(fname,"r")
    totlength = ds.dim["time"] :: Int

    newfnames = String[]

    for (range,name)  in splitindices(totlength,split)
        newfname = replace(fname,".nc" => "." * name * ".nc")
        push!(newfnames,newfname)
        println("$name: indices: $range; fraction: $(length(range)/totlength), $newfname ")
        DINCAE_utils.ncsplit(fname,newfname,time = range)
    end

    return newfnames
end

"""
    mask, count_present = DINCAE_utils.compute_mask(X; minseafrac = 0.05)

Compute the land-sea mask based on the data in `X` (with dimensions
longitude, latitude and time). A grid point is considered land
(i.e. coresponding mask is false) if the fraction of missing data exceeds
`minseafrac`.
"""
function compute_mask(DT; minseafrac = 0.05)
    len = size(DT)[3]
    data = DT[:,:,1]
    count_present = Int.(.!ismissing.(data))

    for i = 2:len
        data = DT[:,:,i]
        count_present += Int.(.!ismissing.(data))
    end

    mask = count_present ./ len .> minseafrac
    return mask,count_present
end


"""
    DINCAE_utils.add_mask(fname,varname; minseafrac = 0.05)

Adds a mask based on minimum fraction of sea data in the NetCDF variable
`varname`.
"""
function add_mask(fname,varname; minseafrac = 0.05)
    Dataset(fname,"a") do ds
        if haskey(ds,"mask")
            @info("mask already present in $fname")
            return
        end

        var = ds[varname]
        mask, count_nomissing = compute_mask(var; minseafrac = minseafrac)
        defVar(ds,"mask",Int8.(mask[:,:,1]),("lon","lat"), attrib = [
            "long_name" => "mask (sea=1, land=0)"
        ])
        defVar(ds,"count_nomissing",Int32.(count_nomissing[:,:,1]),("lon","lat"), attrib = [
            "long_name" => "number of persent data"
        ])
    end
end


"""
    DINCAE_utils.addcvpoint(fname,varname; mincvfrac = 0.10)

Add cross-validation points to a dataset. This functions will withheld data in
the time instance with the highest data coverage (i.e. the cleanest image) 
using the cloud mask from another time instance (choosen at random). Then 
the algorithm will proceed the next image with the second highest coverage and so on.
The algorithm stops when the fraction of cross-validation data point is at minimum 
`mincvfrac` (per default 0.1 or 10%).
"""
function addcvpoint(fname,varname; mincvfrac = 0.10)
    fname_cv = replace(fname,r".nc$" => "_add_clouds.nc")
    cp(fname,fname_cv,force=true)
    n_cv = Int[]

    Dataset(fname_cv,"a") do ds
        data = ds[varname][:,:,:];
        time = ds["time"][:];

        if !haskey(ds,"mask")
            @info "no mask in $fname"
            mask = trues(size(data)[1:2]...,1)
        else
            mask = ds["mask"][:,:][:,:,1:1] .== 1
        end

        data[.!ismissing.(data) .& .!mask] .= missing

        nvalid = sum(.!ismissing.(data))

        ncv = 0

        tmp = data[:,:,1]
        nmissing = sum(ismissing,data,dims=[1,2])[1,1,:]

        for n_dest in sortperm(nmissing)
            n_source = rand(1:size(data,3))

            tmp = data[:,:,n_dest]
            nmissing_before = sum(ismissing,tmp)

            tmp[ismissing.(data[:,:,n_source])] .= missing;
            nmissing_after = sum(ismissing,tmp)

            data[:,:,n_dest] = tmp

            push!(n_cv,n_dest)

            ncv += nmissing_after - nmissing_before

            if ncv >= mincvfrac * nvalid
                break
            end
            #@show n_dest,time[n_dest],nmissing_after - nmissing_before,ncv,mincvfrac * nvalid
        end

        @info("number cross-validation points ",ncv)
        @info("percentage of cross-validation points ",100*ncv/nvalid)

        ds[varname][:] = data
    end
    return fname_cv
end

function listcvimages(case)
    ds_cv = Dataset(case.fname_cv);
    ds_orig = Dataset(case.fname_orig);
    mask = ds_cv["mask"][:,:] .== 1

    image_index = Int[]
    for n = 1:length(ds_cv["time"])
        data_cv = ds_cv[case.varname][:,:,n]
        data_orig = ds_orig[case.varname][:,:,n]
        data_orig[.!mask] .= missing

        ncv = sum(ismissing.(data_cv) .&  .!ismissing.(data_orig))
        if ncv > 0
            push!(image_index,n)
        end
    end

    close(ds_cv)
    close(ds_orig)
    return image_index
end



function seasonalaverage(SST,time; DT = 30, cycle_len = 365)
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

                if count > 0
                    mSST[i,j,nn] /= count
                else
                    mSST[i,j,nn] = missing
                end

#                #for i = 1:2
#                sel .= abs.(  mod.( doy  .- doy[nn] .+ 365/2, 365) .- 365/2) .<= DT÷2;
#                mSST[i,j,nn] = mean(@view SST[i,j,sel])
            end
        end
    end
    return mSST

end


function remove_seasonal_cycle(SST,SSTtime; DT = 30, cycle_len = 365)
    doy = Dates.dayofyear.(SSTtime);
    mSST2 = seasonalaverage(
        SST,SSTtime;
        DT = DT, cycle_len = cycle_len);

    SSTa = similar(SST);
    for n = 1:size(SST,3)
        SSTa[:,:,n] = SST[:,:,n] - mSST2[:,:,doy[n]]
    end
    return SSTa,mSST2
end

function std_around_seasonalaverage(SST,SSTtime; DT = 30, cycle_len = 365)
    doy = Dates.dayofyear.(SSTtime);
    mSST2 = seasonalaverage(
        SST,SSTtime;
        DT = DT, cycle_len = cycle_len);

    SSTa = similar(SST);
    for n = 1:size(SST,3)
        SSTa[:,:,n] = SST[:,:,n] - mSST2[:,:,doy[n]]
    end

    count = sum(.!ismissing.(SSTa), dims = 3)
    SSTa[ismissing.(SST)] .= 0

    SST_std = sqrt.(sum(SSTa.^2,dims = 3) ./ count)[:,:,1];
    return SST_std
end
