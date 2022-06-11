#################### Self Adjusting Feature Detection script - With Indigo parameters ####################

using Distributed
addprocs(8)
# @everywhere push!(LOAD_PATH,$"/Users/saersamanipour/Desktop/dev/ongoing/juliahrms")
#@everywhere push!(LOAD_PATH,pwd())
#@everywhere include("FeatureDetection.jl")
@everywhere using SAFD
#@everywhere using FeatureDetection
@everywhere using MS_Import
@everywhere using ProgressMeter
@everywhere using BenchmarkTools


############################################################################
function file_scanner(path2files)
    Dname = readdir(path2files)
    ind=zeros(length(Dname),1)
    @showprogress "Files are being scanned..." for i=1:length(Dname)
        sleep(0.1)
        if Dname[i][1] != "." && string(Dname[i][1]) != "." && length(Dname[i]) > 5
            if Dname[i][end-4:end] == "mzXML" || Dname[i][end-4:end]== "mzxml" || Dname[i][end-4:end] == "MZXML"
                ind[i]=1
            elseif Dname[i][end-3:end] == "CDF" || Dname[i][end-3:end]== "cdf"
                ind[i]=1
            end
        end
    end
    return(Dname,ind)
end

###########################
function file_package(Dname,ind,n)
    n_pac = ceil(sum(ind)/n)
    Dname_s = Dname[findall(x -> x > 0,ind)]
    filenames_pac =  Array{Any}(undef,Int(n_pac))
    ii = Int.(1:n:n*n_pac+n)
    for i=1:Int(n_pac)
        if ii[i+1]-1 < length(Dname_s)
            filenames_pac[i] = Dname_s[ii[i]:ii[i+1]-1]
        else
            filenames_pac[i] = Dname_s[ii[i]:end]
        end
    end
    return filenames_pac
end

###########################
function par_make(path2files,filenames_pac_s,mz_thresh,Int_thresh)
    mz_vales_mat=Array{Any}(undef,length(filenames_pac_s))
    mz_int_mat=Array{Any}(undef,length(filenames_pac_s))
    t0_mat=Array{Float64}(undef,length(filenames_pac_s))
    t_end_mat=Array{Float64}(undef,length(filenames_pac_s))
    m_mat=Array{Any}(undef,length(filenames_pac_s))
    Rt_mat=Array{Any}(undef,length(filenames_pac_s))
    pathin_mat=Array{String}(undef,length(filenames_pac_s))
    @showprogress "File packages are being created..." for i=1:length(filenames_pac_s)
        sleep(0.1)
        filenames=[filenames_pac_s[i]]
        mz_vals,mz_int,t0,t_end,m,pathin,msModel,msIonisation,msManufacturer,polarity,Rt,centroid=import_files_MS1(path2files,filenames,mz_thresh,Int_thresh)
        mz_vales_mat[i]=mz_vals
        mz_int_mat[i]=mz_int
        t0_mat[i]=t0
        t_end_mat[i]=t_end
        m_mat[i]=m[end-1]
        pathin_mat[i]=path2files
        Rt_mat[i]=Rt
    end
    return (mz_vales_mat,mz_int_mat,t0_mat,t_end_mat,m_mat,pathin_mat,Rt_mat)
end


###########################
function par_make_cent(path2files,filenames_pac_s,mz_thresh,Int_thresh)
    mz_vales_mat=Array{Any}(undef,length(filenames_pac_s))
    mz_int_mat=Array{Any}(undef,length(filenames_pac_s))
    t0_mat=Array{Float64}(undef,length(filenames_pac_s))
    t_end_mat=Array{Float64}(undef,length(filenames_pac_s))
    m_mat=Array{Any}(undef,length(filenames_pac_s))
    Rt_mat=Array{Any}(undef,length(filenames_pac_s))
    pathin_mat=Array{String}(undef,length(filenames_pac_s))
    #dm_c_mat = Array{String}(undef,length(filenames_pac_s))
    @showprogress "File packages are being created..." for i=1:length(filenames_pac_s)
        sleep(0.1)
        filenames=[filenames_pac_s[i]]
        #chrom=import_files(path2files,filenames,mz_thresh,Int_thresh)
        mz_vals,mz_int,t0,t_end,m,pathin,msModel,msIonisation,msManufacturer,polarity,Rt,centroid=import_files_MS1(path2files,filenames,mz_thresh,Int_thresh)
        #mz_val_cent,mz_int_cent,dm_c = centroid(chrom["MS1"]["Mz_values"],chrom["MS1"]["Mz_intensity"],min_int,res);
        mz_vales_mat[i] = mz_vals
        mz_int_mat[i] = mz_int
        t0_mat[i]=t0
        t_end_mat[i]=t_end
        m_mat[i]=m[end-1]
        pathin_mat[i]=path2files
        Rt_mat[i]=Rt
        #dm_c_mat[i] = dm_c
    end
    return (mz_vales_mat,mz_int_mat,t0_mat,t_end_mat,m_mat,pathin_mat,Rt_mat)
end


###########################
function pmap_centroid(mz_vales_mat,mz_int_mat,res,min_int)
    Min_int=repeat([min_int],size(mz_vales_mat,1))
    Res=repeat([res],size(mz_vales_mat,1))
    mz_vales_mat_cent = Array{Any}(undef,size(mz_vales_mat,1))
    mz_int_mat_cent = Array{Any}(undef,size(mz_vales_mat,1))
    cdm_mat = Array{Any}(undef,size(mz_vales_mat,1))
    output = pmap(centroid,mz_vales_mat,mz_int_mat,Min_int,Res;retry_delays = zeros(2))
    GC.gc()
    for i =1:size(mz_vales_mat,1)
        mz_vales_mat_cent[i] = output[i][1]
        mz_int_mat_cent[i] = output[i][2]
        cdm_mat[i] = output[i][3]
    end
    return(mz_vales_mat_cent,mz_int_mat_cent,cdm_mat)
end


###########################
function pmap_processing(mz_vales_mat,mz_int_mat,m_mat,pathin_mat,Rt_mat,
    max_numb_iter,max_t_peak_w,res,min_ms_w,r_thresh,min_int,sig_inc_thresh,S2N,min_peak_w_s)
    Max_it=repeat([max_numb_iter],size(pathin_mat,1))
    Max_t_peak_w=repeat([max_t_peak_w],size(pathin_mat,1))
    Res=repeat([res],size(pathin_mat,1))
    Min_ms_w=repeat([min_ms_w],size(pathin_mat,1))
    R_thresh=repeat([r_thresh],size(pathin_mat,1))
    Min_int=repeat([min_int],size(pathin_mat,1))
    Sig_inc_thresh=repeat([sig_inc_thresh],size(pathin_mat,1))
    S2N1=repeat([S2N],size(pathin_mat,1))
    Min_peak_w_s=repeat([min_peak_w_s],size(pathin_mat,1))
    pmap(safd_s3D,mz_vales_mat,
    mz_int_mat,Rt_mat,m_mat,pathin_mat,Max_it,Max_t_peak_w,Res,
    Min_ms_w,R_thresh,Min_int,Sig_inc_thresh,S2N1,Min_peak_w_s;retry_delays = zeros(2))
    GC.gc()
end


###########################
function pmap_processing_cent(mz_vales_mat,mz_int_mat,m_mat,pathin_mat,Rt_mat,dm_c_mat,
    max_numb_iter,max_t_peak_w,res,min_ms_w,r_thresh,min_int,sig_inc_thresh,S2N,min_peak_w_s,method,dm)
    Max_it=repeat([max_numb_iter],size(pathin_mat,1))
    Max_t_peak_w=repeat([max_t_peak_w],size(pathin_mat,1))
    Res=repeat([res],size(pathin_mat,1))
    Min_ms_w=repeat([min_ms_w],size(pathin_mat,1))
    R_thresh=repeat([r_thresh],size(pathin_mat,1))
    Min_int=repeat([min_int],size(pathin_mat,1))
    Sig_inc_thresh=repeat([sig_inc_thresh],size(pathin_mat,1))
    S2N1=repeat([S2N],size(pathin_mat,1))
    Min_peak_w_s=repeat([min_peak_w_s],size(pathin_mat,1))
    method_mat = repeat([method],size(pathin_mat,1))
    dm_mat = repeat([dm],size(pathin_mat,1))
    if method == "BG"
        pmap(safd_s3d_cent,mz_vales_mat,mz_int_mat,Rt_mat,m_mat,pathin_mat,Max_it,
        Max_t_peak_w,Res,Min_ms_w,R_thresh,Min_int,Sig_inc_thresh,S2N1,Min_peak_w_s,method_mat;retry_delays = zeros(2))
        GC.gc()
    elseif method == "CT"
        pmap(safd_s3d_cent,mz_vales_mat,mz_int_mat,Rt_mat,m_mat,pathin_mat,Max_it,
        Max_t_peak_w,Res,Min_ms_w,R_thresh,Min_int,Sig_inc_thresh,S2N1,Min_peak_w_s,method_mat,dm_mat;retry_delays = zeros(2))
        GC.gc()
    elseif method == "RFM"
        pmap(safd_s3d_cent,mz_vales_mat,mz_int_mat,Rt_mat,m_mat,pathin_mat,Max_it,
        Max_t_peak_w,Res,Min_ms_w,R_thresh,Min_int,Sig_inc_thresh,S2N1,Min_peak_w_s,method_mat;retry_delays = zeros(2))
        GC.gc()
    elseif method == "MDM"
        pmap(safd_s3d_cent,mz_vales_mat,mz_int_mat,Rt_mat,m_mat,pathin_mat,Max_it,
        Max_t_peak_w,Res,Min_ms_w,R_thresh,Min_int,Sig_inc_thresh,S2N1,Min_peak_w_s,method_mat,dm_c_mat;retry_delays = zeros(2))
        GC.gc()
    end
end


###########################
function safd_pmap(path2files,n,Int_thresh,mz_thresh,max_numb_iter,
    max_t_peak_w,res,min_ms_w,r_thresh,min_int,sig_inc_thresh,S2N,min_peak_w_s)
    Dname,ind = file_scanner(path2files)
    filenames_pac = file_package(Dname,ind,n)
    @showprogress "File packages are being proccessed..." for i=1:size(filenames_pac,1)
        sleep(0.1)
        mz_vales_mat,mz_int_mat,t0_mat,t_end_mat,m_mat,pathin_mat,Rt_mat = par_make(path2files,filenames_pac[i],mz_thresh,Int_thresh)
        println("The file package is imported and will be proccessed.")
        pmap_processing(mz_vales_mat,mz_int_mat,m_mat,pathin_mat,Rt_mat,
        max_numb_iter,max_t_peak_w,res,min_ms_w,r_thresh,min_int,sig_inc_thresh,S2N,min_peak_w_s)
        mz_vales_mat = []
        mz_int_mat = []
        t0_mat = []
        t_end_mat = []
        m_mat = []
        pathin_mat = []
        Rt_mat = []
        GC.gc()
    end
end


###########################
function safd_pmap_cent(path2files,n,Int_thresh,mz_thresh,max_numb_iter,
    max_t_peak_w,res,min_ms_w,r_thresh,min_int,sig_inc_thresh,S2N,min_peak_w_s,method,dm)
    Dname,ind = file_scanner(path2files)
    filenames_pac = file_package(Dname,ind,n)
    @showprogress "File packages are being proccessed..." for i=1:size(filenames_pac,1)
        sleep(0.1)
        mz_vales_mat,mz_int_mat,t0_mat,t_end_mat,m_mat,pathin_mat,Rt_mat = par_make_cent(path2files,filenames_pac[i],mz_thresh,Int_thresh)
        println("The file package is being converted to centroids.")
        mz_vales_mat_cent,mz_int_mat_cent,cdm_mat = pmap_centroid(mz_vales_mat,mz_int_mat,res,min_int)
        pmap_processing_cent(mz_vales_mat_cent,mz_int_mat_cent,m_mat,pathin_mat,Rt_mat,cdm_mat,
        max_numb_iter,max_t_peak_w,res,min_ms_w,r_thresh,min_int,sig_inc_thresh,S2N,min_peak_w_s,method,dm)
        mz_vales_mat = []
        mz_int_mat = []
        t0_mat = []
        t_end_mat = []
        m_mat = []
        pathin_mat = []
        Rt_mat = []
        GC.gc()
    end
end


##############################################################################
n=1


path2files = "/media/saer/Elements SE/Data Jiani/"
#filenames=["filename.mzXML"]
mz_thresh = [0,600]
Int_thresh = 8000

# Feature detection parameters
max_numb_iter = 5000
max_t_peak_w = 300 # or 30
res = 30000
min_ms_w = 0.02
r_thresh = 0.75
min_int = 10000
sig_inc_thresh = 5
S2N = 2
min_peak_w_s = 3

dm = 10
method = "MDM" # "BG" # "RFM" # "CT" #

# using SAFD
safd_pmap(path2files,n,Int_thresh,mz_thresh,max_numb_iter,
   max_t_peak_w,res,min_ms_w,r_thresh,min_int,sig_inc_thresh,S2N,min_peak_w_s)

#@time safd_pmap_cent(path2files,n,Int_thresh,mz_thresh,max_numb_iter,
#    max_t_peak_w,res,min_ms_w,r_thresh,min_int,sig_inc_thresh,S2N,min_peak_w_s,method,dm)

folder = "Full/SS3 +/"
files = readdir("/media/saer/Elements SE/Data Jiani/"*folder)

for f in files
    mv("/media/saer/Elements SE/Data Jiani/"*folder*f,
       "/media/saer/Elements SE/Data Jiani/"*folder*"SS3_pos_"*f)
end

for f in files
    p = split(f, ".")
    if length(p) .== 2
        # println(p)
        continue
    elseif length(p) .> 3
        println("found more dots")
        continue
    else
        # println(" fixing  " )
    end
    println(f)
    mv("/media/saer/Elements SE/Data Jiani/"*folder*f,
       "/media/saer/Elements SE/Data Jiani/"*folder*p[1]*"-"*p[2]*"."*p[3])
end

#######################################################################

using SAFD
using CSV
using DataFrames

#path2files_a = "/media/saer/Elements SE/Data Jiani/SS1 pos full time" # change the path to the locat of the files that you want to align
#path2files_a = "/media/saer/Elements SE/Data Jiani/SS2 pos full time"
#path2files_a = "/media/saer/Elements SE/Data Jiani/SS3 pos full time"

#path2files_a = "/media/saer/Elements SE/Data Jiani/SS1 neg full time" # change the path to the locat of the files that you want to align
#path2files_a = "/media/saer/Elements SE/Data Jiani/SS2 neg full time"
#path2files_a = "/media/saer/Elements SE/Data Jiani/SS3 neg full time"

path2files_a = "/media/saer/Elements SE/Data Jiani/SS_complete_neg_full_time"
#path2files_a = "/media/saer/Elements SE/Data Jiani/SS_complete_pos_full_time"

int_table, a_table, r_table = feature_align(path2files_a)
CSV.write("/media/saer/Elements SE/Data Jiani/SS_complete_neg_full_time.csv",int_table)
println("end")
