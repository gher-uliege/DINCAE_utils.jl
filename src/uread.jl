
function uread(filename)
    f = open(filename)

    # ignore header
    header = zeros(Int32,20)
    read!(f,header)

    reclen1 = read(f,Int32)

    if reclen1 == 24
        # no swapping is necessary
        doswap = identity
    else
        reclen1 = bswap(reclen1)
        doswap = bswap
    end
    @assert reclen1 == 24

    imax = doswap(read(f,Int32))
    jmax = doswap(read(f,Int32))
    kmax = doswap(read(f,Int32))

    iprec = doswap(read(f,Int32))
    nbmots = doswap(read(f,Int32))
    valex = doswap(read(f,Float32))

    reclen2 = doswap(read(f,Int32))

    @assert reclen1 == reclen2

    nl = (imax*jmax*kmax) ÷ nbmots
    ir = imax*jmax*kmax - nbmots * nl;

    ide = 1;

    T = (iprec == 4 ? Float32 : Float64 )

    data = zeros(T,imax,jmax,kmax)

    tmp = zeros(T,nbmots)

    for i = 1:nl
        #read(f,@view data[ide:ide+nbmots])
        reclen1 = doswap(read(f,Int32))
        read!(f,tmp)
        data[ide:ide+nbmots-1] = doswap.(tmp)
        ide = ide + nbmots
        reclen2 = doswap(read(f,Int32))
        @assert reclen1 == reclen2
    end

    tmp = zeros(T,ir)

    reclen1 = doswap(read(f,Int32))
    read!(f,tmp)
    data[ide:ide+ir-1] = doswap.(tmp)
    reclen2 = doswap(read(f,Int32))
    @assert reclen1 == reclen2

    return data,T(valex)
end
#read(f,@view data[ide:ide+ir-1])

function gread(filename,val)
    data,valex = uread(filename)
    return replace(data,valex => val)
end
