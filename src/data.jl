"""
    DINCAE_utils.splitdata(fname,split = [("train",0.7),("dev",0.2),("test",0.1)])

Split the NetCDF file `fname` into training, developpement and test dataset
with the specified fraction.
"""
function splitdata(fname,split = [("train",0.7),("dev",0.2),("test",0.1)])
    ds = NCDataset(fname,"r")
    totlength = ds.dim["time"] :: Int

    i = 0
    newfnames = String[]

    for (name,fraction)  in split
        len = round(Int,totlength * fraction)
        range =  (i+1) : min(i+len, totlength)
        i += len
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
the time slices with the highest data coverage using the data mask from other
time instances.

Adds a mask based on minimum fraction of sea data in the NetCDF variable
`varname`.
"""
function addcvpoint(fname,varname; mincvfrac = 0.10)
    fname_cv = replace(fname,r".nc$" => "_add_clouds.nc")
    cp(fname,fname_cv,force=true)
    n_cv = Int[]

    Dataset(fname_cv,"a") do ds
        data = ds[varname][:,:,:];
        time = ds["time"][:];
        mask = ds["mask"][:,:][:,:,1:1] .== 1

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
            @show n_dest,time[n_dest],nmissing_after - nmissing_before,ncv,mincvfrac * nvalid
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
    for n = 1:ds_cv.dim["time"]
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

