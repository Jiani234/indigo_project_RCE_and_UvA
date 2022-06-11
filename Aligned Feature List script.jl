#################### Aligned Feature List script ####################  

# __precompile__()

# module FeatureAlign

#push!(LOAD_PATH,pwd())
using Base

using Plots
using CSV
using DataFrames
using Glob
using Statistics
using XLSX

# export feature_alignment_wraper

#####
# Starting the Functions

###
# A function to import the files

function report_import(path2files,ext,source)
    if source == "Internal"
        nm=glob("*_report_*."*ext,path2files)
        Reps=Dict()
        name=[]
        for i=1:size(nm,1)
            nf=splitdir(nm[i])
            nf1=split(nf[2],".")
            name=vcat(name,nf1[1])
            df = CSV.File(nm[i]) |> DataFrame
            Reps[nf1[1]]=df
            df=[]
        end
    end

    return(Reps,name)

end

# Reps,name=report_import(path2files,ext,source)

###
# Feature Creation

function feature_creation(Reps,name)

    features = Reps[name[1]]

    for i =2:length(name)
        features=vcat(features,Reps[name[i]])

    end

    sort!(features,[:MeasMass,:ScanNum])

    return features

end

# feature_list=feature_creation(Reps,name)

###
# Master feature list creation

function master_feature_list(features)

    m_feature_list = zeros(size(features,1),9)

    features_temp = deepcopy(features)
    features_temp.Nr = 1:size(features,1)

    for i =1:size(features,1)
        if features_temp.Nr[i] == 0
            continue
        end
        println(i)
        win = [features_temp.ScanNum[i]- 2*features_temp.ScanInPeak[i],features_temp.ScanNum[i]+ 2*features_temp.ScanInPeak[i]]
        if win[1] <= 0
            win[1] =1
        elseif win[2] > maximum(features_temp.ScanNum)
            win[2] = maximum(features_temp.ScanNum)
        end
        p_feat = features_temp[findall(x -> win[1] <= x <= win[2],features_temp.ScanNum),:]
        m_dev = abs.(p_feat.MeasMass[1] .- p_feat.MeasMass)

        p_feat_1 = p_feat[m_dev .< 1,:]
        m_tol = round(maximum(p_feat_1.MaxMass .- p_feat_1.MinMass),digits=3)
        m_e = abs.(p_feat_1.MeasMass[1] .- p_feat_1.MeasMass)
        s_feat = p_feat_1[m_e .<= m_tol,:]

        m_feature_list[i,1] = minimum(s_feat.ScanNum)
        m_feature_list[i,2] = maximum(s_feat.ScanNum)
        m_feature_list[i,3] = mean(s_feat.ScanNum)

        m_feature_list[i,4] = minimum(s_feat.RtStart)
        m_feature_list[i,5] = maximum(s_feat.RtEnd)
        m_feature_list[i,6] = mean(s_feat.Rt)

        m_feature_list[i,7] = minimum(s_feat.MinMass)
        m_feature_list[i,8] = maximum(s_feat.MaxMass)
        m_feature_list[i,9] = mean(s_feat.MeasMass)

        features_temp[s_feat.Nr,:] .= 0

    end

    m_feature_list_f = m_feature_list[m_feature_list[:,1] .>0,:]


    return m_feature_list_f

end

# m_feature_list=master_feature_list(feature_list)

##########################
# Feature find
# i=3

function feature_find(file,m_feature_list_f)
    int=zeros(size(m_feature_list_f,1),1)
    # area=zeros(size(m_feature_list_f,1),1)
    # res=zeros(size(m_feature_list_f,1),1)

    for i=1:size(m_feature_list_f,1)
        tv1 = m_feature_list_f[i,:]

        ind = intersect(findall(x -> tv1[2] >= x >= tv1[1],file.ScanNum), findall(x -> tv1[8] >= x >= tv1[7],file.MeasMass))

        if length(ind)>1
            #println(i)

            int[i]=maximum(file.Int[ind])
            # area[i]=sum(file.Area[ind])
            # res[i]=mean(file.MediRes[ind])

        elseif length(ind)==1
            #println(i)

            int[i]=file.Int[ind][1]
            # area[i]=file.Area[ind][1]
            # res[i]=file.MediRes[ind][1]
        end

    end

    return int#,area,res)

end

# int,area,res=feature_find(Reps,Name,m_feature_list)

###
# Alignment function

function feature_alignment(Reps,name,m_feature_list_f)

    Inten=zeros(size(m_feature_list_f,1),size(name,1));
    # Area=zeros(size(m_feature_list_f,1),size(name,1));
    # Ress=zeros(size(m_feature_list_f,1),size(name,1));

    for i=1:size(name,1)
        println(i)
        file=Reps[name[i]]
        # int,area,res=feature_find(file,m_feature_list_f)
        int=feature_find(file,m_feature_list_f)

        #scatter(area,int)


        Inten[:,i]=int
        # Area[:,i]=area
        # Ress[:,i]=res
        #scatter(Area[:,i],Inten[:,i])

    end

    table1=DataFrame(hcat(1:size(m_feature_list_f,1),m_feature_list_f),[:Nr,:MinScanNum,:MaxScanNum,:AveScanNum,:MinRt,:MaxRt,:AveRt,
    :MinMass,:MaxMass,:AveMass])
    table1.Nr = 1:size(m_feature_list_f,1)

    table_int=DataFrame(Inten, :auto)
    # table_are=DataFrame(Area, :auto)
    # table_res=DataFrame(Ress, :auto)

    rename!(table_int,Symbol.(name))
    # rename!(table_are,Symbol.(name))
    # rename!(table_res,Symbol.(name))

    table_int_f=hcat(table1,table_int)
    # table_are_f=hcat(table1,table_are)
    # table_res_f=hcat(table1,table_res)

    sort!(table_int_f,[:AveScanNum,:AveMass])
    # sort!(table_are_f,[:AveScanNum,:AveMass])
    # sort!(table_res_f,[:AveScanNum,:AveMass])

    return table_int_f#,table_are_f,table_res_f)

end

# table_int_f,table_are_f,table_res_f=feature_alignment(Reps,name,m_feature_list)
#########################################
function feature_alignment_wraper(path2files,ext,source)

    Reps,name=report_import(path2files,ext,source)

    println("All the files have been imported.")

    feature_list=feature_creation(Reps,name)

    println("All possible features are created.")

    m_feature_list_f=master_feature_list(feature_list)
    m=size(m_feature_list_f,1)
    #m=1000
    m1=string(m)
    println("A master feature list of "*m1 * " features has been created.")
    #int,area,res=feature_find(Reps,Name,m_feature_list)
    # table_int_f,table_are_f,table_res_f=feature_alignment(Reps,name,m_feature_list_f)
    table_int_f=feature_alignment(Reps,name,m_feature_list_f)

    output=joinpath(path2files,"FeatureList_Aligned.xlsx")

    range_i = Int(round(size(table_int_f,1)/3))
    CSV.write(joinpath(path2files,"FeatureList_Aligned_Intensities_1.csv"),table_int_f[1:range_i,:])
    CSV.write(joinpath(path2files,"FeatureList_Aligned_Intensities_2.csv"),table_int_f[range_i+1:2*range_i,:])
    CSV.write(joinpath(path2files,"FeatureList_Aligned_Intensities_3.csv"),table_int_f[2*range_i+1:end,:])
    # CSV.write(joinpath(path2files,"FeatureList_Aligned_Areas.csv"),table_are_f)
    # CSV.write(joinpath(path2files,"FeatureList_Aligned_Resolutions.csv"),table_res_f)

    # XLSX.writetable(output, Intensities=(collect(DataFrames.eachcol(table_int_f)),
    #  DataFrames.names(table_int_f)  ), Areas=(collect(DataFrames.eachcol(table_are_f)),
    #   DataFrames.names(table_are_f)  ), Resolutions=(collect(DataFrames.eachcol(table_res_f)),
    #    DataFrames.names(table_res_f)  ), overwrite=true)


    println("The final report has been saved in the output path")

    return table_int_f

end

# end #End of the module
using SAFD

path = "C:\\Users\\dherwer\\OneDrive - UvA\\Projects\\Jiani\\Feature lists\\Final major feature alignment"
ACTUALpath2files = "/media/saer/Elements SE/Data Jiani/pos full time/"

files = readdir(path)
for f in files
    path2files = path*"\\"*f
    #path2files = "/Volumes/SAER HD/Data/Temp_files/Phil/Test_align/"
    mod = readdir(path2files)
    for k in mod
        ACTUALpath2files = path2files*"\\"*k
        v_n_mz = "MeasMass"
        v_n_rt = "Rt"
        v_n_int = "Int"
        mz_tol = 0.02
        rt_tol = 0.1

        # Internal
        feature_align(ACTUALpath2files)
    end
end
# Internal with tols
# feature_align(path2files,mz_tol,rt_tol)
# External
# feature_align(path2files,mz_tol,rt_tol,v_n_mz,v_n_rt,v_n_int)

####

path2files="/Users/jianihu/Documents/Sample set 3 Compleet/Sample set 3 na SAFD Asb 24-80 negatief en positief/Positief AsB 24-80"
ext="csv"
source="Internal"

@time feature_alignment_warper(path2files,ext,source)

path2files="/Users/jianihu/Documents/Sample set 3 Compleet/Sample set 3 na SAFD Asb 24-80 negatief en positief/Positief AsB 24-80"
ext="csv"
source="Internal"
