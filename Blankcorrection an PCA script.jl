##### Blank correction & PCA SVD and Eigen Decomposition script #####

## needed packages
using DataFrames, XLSX, Statistics, Distributions, CSV, LinearAlgebra, Plots

## import FeatureList
ssc_neg_Intensities = DataFrame(XLSX.readtable("/Users/jianihu/Documents/FeatureList_Aligned_Full_Neg.xlsx","Intensities")...)      # Data import
ssc_neg_Intensities = ssc_neg_Intensities[ 3 .< ssc_neg_Intensities[!,"AverageRT"] .<30 ,:]                                         # Selecting chromatogram data between 3.0 - 30 min

## file in sample information 
Sample = "SS All"       # Sample set 
Charge = "Negative"     # Negative or Positive mode 
plotv = "V1"            # version of plots 
T1 = "PCA"              # PCA or PLSDA

## Matrix clean up
# Blanco signal correction
function blanco_clean_up(data_mat)                      # data_mat = ss1_pos_Intensities of een andere
    blanco = contains.(names(data_mat),"Blanco")        # Blanco eruit filteren voor de correctie
    standaard = contains.(names(data_mat),"Coch")       # Coch calibration sample filter
    ADD = contains.(names(data_mat),"ADD")              # ADD sample filter
    MME = contains.(names(data_mat),"MME")
    sLPE = contains.(names(data_mat),"sLPE")
    SS1 = contains.(names(data_mat),"SS1")
    SS2 = contains.(names(data_mat),"SS2")
    SS3 = contains.(names(data_mat),"SS3")
    Margriet = contains.(names(data_mat),"M ")
    ACIAR = contains.(names(data_mat),"ACIAR")
    WrB = contains.(names(data_mat),"WrB")
    trade = contains.(names(data_mat),"trade")
    
    samples = contains.(names(data_mat),"SS") .& (.! blanco) .& (.! standaard) # First part selecting the data and second part erasing/converting the other data to 0 
    
    m = mean(Matrix(data_mat[:, blanco]), dims=2)           # Mean berekening van de 2e dims --> alle rijen

    for i = 1:size(data_mat,2)                              # Telt het aantal columns van data_mat, met 2 geven we de dims aan (voor size is andere manier om dims aan te geven)
        if samples[i] == 1                                  # Gaat kijken of er een samples aanwezig is => 1 = true en 0 = false
            data_mat[vec(data_mat[:, i] .< 3*m), i].= 0     # Kijken door de hele rij van de bijbehorende i of het getal < 3xm is, zo ja dan wordt dit een 0, want die nemen we niet mee
        end                                                 # want het signaal moet minimaal 3x de blanco zijn
    end

    data_mat_b_clean = data_mat[:,sortperm(names(data_mat))]        # Sorteert de namen op abc. van de df, nu nummers ook bijelkaar
    naam = names(data_mat_b_clean)                                  # Vector met alle namen) 

    v = sum(Matrix(data_mat[!,samples[1:end]]).> 0, dims=2)         # alle samlples bekijken > 0 --> naar 0

    for u = 1:size(data_mat_b_clean,2)
        if (samples[u] == 1)
            data_mat[!,u][vec(v .< length(samples)*0.10)].= 0       # <2 = less than 50% present --> 0
        end
    end

    for i = 1:size(data_mat_b_clean,2)                          # Zelfde als vorige For statement
        if (samples[i] == 1) .& (!contains.(naam[i],"Coch"))    # Coch standaarden eruit filteren
            println(naam[i])                                                
            sn = split(split(naam[i]," ")[1],"_")[end]          # sn = sample naam, AsB gedeelte selecteren --> split de naam bij de spatie en selecteerd het 1e stuk en daarna op _ filteren en het laatste stuk (AsB) slecteren
            snum = split(split(naam[i]," ")[2],"-")[1]          # snum = sample nummer selecteren, zelfde principe
            sel = []                                            # Lege vector creeeren
            tempnm = naam[findall(contains.(naam,sn))]          # sample type filteren

            for k in tempnm
                if split(k[length(split(naam[i]," ")[1])+2:end],"-")[1] == snum     # Het nummer selecteren van de sample
                    sel =[sel;k]                                                    # sel vullen met sel ; k, dus vector met alle samples die hetzelfde zijn
                end
            end

            t = sum(Matrix(data_mat[!,sel[1:end]]).> 0, dims=2)     # geselecteerde kolomen t= true false van elke kolom/mz waarde

            for s in sel
                data_mat[!,s][vec(t .< 2)].= 0                      # <2 = less than 50% present --> 0
            end
        end
    end

    return data_mat, samples

end

z,samples = blanco_clean_up(ssc_neg_Intensities)        # z,samples zorgt ervoor dat de  
z = z[vec(any(Matrix(z)[:,samples].!=0,dims=2)),:]      # alle rijen, die gelijk zijn aan -0 weg

## Save import data after blank correction 
SSC_neg = CSV.write("/Users/jianihu/Documents/PCA full negative/FeatureList_Aligned_$Charge _full_PLSDA $Sample $plotv.csv", z[:,samples])

b = Matrix(z[:,samples])'
cov(b)

function standard(b)            # Normalise the data
    m = mean(b,dims=1)          # mean berekenen
    stdev = std(b,dims=1)       # stdev berekenen
    bs = (b .- m) ./ stdev      # van alles de mean aftrekken en delen door de stdev => meancentering & scaling
    return bs
end

## binaire matrix ##                                            
y_pca = zeros(size(b,1),length(unique(y)))
yvalues = unique(y)
for i = 1:length(unique(y))
    y_pca[y.== yvalues[i],i] .= 1
end

b[b .> 0].= 1
                                                
## PCA Eigen Decomposition ##

function pca2(bs)                                   # PCA function --> x = input data en pec = proportion of eigenvallue
    cov_mat=cov(bs)
    value,vec_mat=eigen(cov_mat)                    # calculate eigenvalue & eigen matrix
    indices= sortperm(value ,rev=true)
    explainvariance=(value[indices]) ./sum(value)   # sort based on explained variance
    eigenvecs = vec_mat[:,indices]
    loading = eigenvecs
    scores = (eigenvecs' * bs')'                   
    return explainvariance, loading, scores         # local to global 
end

b = Matrix(z[:,samples])'
bs = standard(b)

explainvariance, loading, scores = pca2(bs)

## explainvariance plot
explainvariance*100
PC1 = round(decimal(explainvariance[1]), digits = 4)*100
PC2 = round(decimal(explainvariance[2]), digits = 4)*100
PC3 = round(decimal(explainvariance[3]), digits = 4)*100
PC4 = round(decimal(explainvariance[4]), digits = 4)*100
                                            
plot_explainV =                                     # plotting explained variace in barplot 
    bar(explainvariance[1:4])
    bar!(frame =:box, dpi = 120, legend = false)
    bar!(title = "$Sample $Charge Explained variance", titlefont = 14)
    bar!(xlabel = "Principle Component", ylabel = "Variance explained (%)", tickfont=font(12),guidefont=font(12))
    savefig("/Users/jianihu/Documents/PCA full negative/Explainedvariance plot $Sample $Charge $plotv.pdf")

## PCA SVD ##  
pca = PCA(n_components=5)
println("start")
pca.fit(b)

variance_e = pca.explained_variance_ratio_

scores = pca.fit_transform(b)
loadings = pca.components_

# explainvariance, loading, scores = pca2(b)
println("end")

pca100 = pca.explained_variance_ratio_*100
PC1 = round(decimal(pca.explained_variance_ratio_[1]), digits = 4)*100
PC2 = round(decimal(pca.explained_variance_ratio_[2]), digits = 4)*100
PC3 = round(decimal(pca.explained_variance_ratio_[3]), digits = 4)*100
PC4 = round(decimal(pca.explained_variance_ratio_[4]), digits = 4)*100

barplot = 
bar(pca100[1:4])
bar!(frame =:box, dpi = 120, legend = false)
bar!(title = "$Sample $Charge Explained variance")
bar!(xlabel = "Principle Component", ylabel = "Variance explained (%)")
savefig("$path2file/Explainedvariance plot $Charge $Sample $plotv.pdf")

## Plotting & Classification ##                                              
meting = names(z)[samples]                          # classifying the different plant species and extraction methods 
y = zeros(size(meting))                             # total 23 variaties
y[contains.(lowercase.(meting),"_asb 1-")] .= 1
y[contains.(lowercase.(meting),"_asb 2-")] .= 1
y[contains.(lowercase.(meting),"_asb 3-")] .= 1
y[contains.(lowercase.(meting),"_asb 4-")] .= 2
y[contains.(lowercase.(meting),"_asb 5-")] .= 2
y[contains.(lowercase.(meting),"_asb 6-")] .= 2
y[contains.(lowercase.(meting),"_asb 7-")] .= 2
y[contains.(lowercase.(meting),"_asb 8-")] .= 2
y[contains.(lowercase.(meting),"_asb 9-")] .= 3
y[contains.(lowercase.(meting),"_asb 10")] .= 3
y[contains.(lowercase.(meting),"_asb 11")] .= 3
y[contains.(lowercase.(meting),"_asb 12")] .= 4
y[contains.(lowercase.(meting),"_asb 13")] .= 4
y[contains.(lowercase.(meting),"_asb 14")] .= 4
y[contains.(meting, "M ")] .= 8
y[contains.(lowercase.(meting),"_asb 15")] .= 5
y[contains.(lowercase.(meting),"_asb 16")] .= 6
y[contains.(lowercase.(meting),"_asb 17")] .= 7
y[contains.(lowercase.(meting),"_asb 18")] .= 7
y[contains.(lowercase.(meting),"_asb 19")] .= 9
y[contains.(lowercase.(meting),"_asb 20")] .= 9
y[contains.(lowercase.(meting),"_asb 21")] .= 9
y[contains.(lowercase.(meting),"_asb 22")] .= 9
y[contains.(lowercase.(meting),"_asb 23")] .= 9
y[contains.(meting, "HallTex")] .= 10
y[contains.(meting, "ACIAR")] .= 11
y[contains.(lowercase.(meting), "_wrb 1-")] .= 12
y[contains.(lowercase.(meting), "_wrb 2-")] .= 12
y[contains.(lowercase.(meting), "_wrb 3-")] .= 12
y[contains.(lowercase.(meting), "_wrb 4")] .= 12
y[contains.(lowercase.(meting), "_wrb 5")] .= 12
y[contains.(lowercase.(meting), "_wrb 6")] .= 12
y[contains.(lowercase.(meting), "_wrb 7")] .= 12
y[contains.(lowercase.(meting), "_wrb 8")] .= 12
y[contains.(lowercase.(meting), "_wrb 9")] .= 12
y[contains.(lowercase.(meting), "_wrb 10")] .= 12
y[contains.(lowercase.(meting), "_wrb 11")] .= 12
y[contains.(lowercase.(meting), "_wrb 12")] .= 12
y[contains.(lowercase.(meting), "_wrb 13")] .= 12
y[contains.(lowercase.(meting), "_wrb 14")] .= 12
y[contains.(lowercase.(meting), "_wrb 15")] .= 12
y[contains.(lowercase.(meting), "_wrb 16")] .= 12
y[contains.(meting, "trade")] .= 13
y[contains.(lowercase.(meting),"_asb 24")] .= 14
y[contains.(lowercase.(meting),"_asb 25")] .= 14
y[contains.(lowercase.(meting),"_asb 26")] .= 14
y[contains.(lowercase.(meting),"_asb 27")] .= 14
y[contains.(lowercase.(meting),"_asb 28")] .= 14
y[contains.(lowercase.(meting),"_asb 29")] .= 14
y[contains.(lowercase.(meting),"_asb 30")] .= 14
y[contains.(lowercase.(meting),"_asb 31")] .= 14
y[contains.(lowercase.(meting),"_asb 32")] .= 14
y[contains.(lowercase.(meting),"_asb 33")] .= 15
y[contains.(lowercase.(meting),"_asb 34")] .= 15
y[contains.(lowercase.(meting),"_asb 35")] .= 15
y[contains.(lowercase.(meting),"_asb 36")] .= 16
y[contains.(lowercase.(meting),"_asb 37")] .= 16
y[contains.(lowercase.(meting),"_asb 38")] .= 16
y[contains.(lowercase.(meting),"_asb 39")] .= 15
y[contains.(lowercase.(meting),"_asb 40")] .= 15
y[contains.(lowercase.(meting),"_asb 41")] .= 15
y[contains.(lowercase.(meting),"_asb 42")] .= 15
y[contains.(lowercase.(meting),"_asb 43")] .= 15
y[contains.(lowercase.(meting),"_asb 44")] .= 15
y[contains.(lowercase.(meting),"_asb 45")] .= 16
y[contains.(lowercase.(meting),"_asb 46")] .= 16
y[contains.(lowercase.(meting),"_asb 47")] .= 16
y[contains.(lowercase.(meting),"_asb 48")] .= 16
y[contains.(lowercase.(meting),"_asb 49")] .= 16
y[contains.(lowercase.(meting),"_asb 50")] .= 16
y[contains.(lowercase.(meting),"_asb 51")] .= 17
y[contains.(lowercase.(meting),"_asb 52")] .= 17
y[contains.(lowercase.(meting),"_asb 53")] .= 17
y[contains.(lowercase.(meting),"_asb 54")] .= 17
y[contains.(lowercase.(meting),"_asb 55")] .= 17
y[contains.(lowercase.(meting),"_asb 56")] .= 17
y[contains.(lowercase.(meting),"_asb 57")] .= 15
y[contains.(lowercase.(meting),"_asb 58")] .= 15
y[contains.(lowercase.(meting),"_asb 59")] .= 15
y[contains.(lowercase.(meting),"_asb 60")] .= 18
y[contains.(lowercase.(meting),"_asb 61")] .= 18
y[contains.(lowercase.(meting),"_asb 62")] .= 18
y[contains.(lowercase.(meting),"_asb 63")] .= 18
y[contains.(lowercase.(meting),"_asb 64")] .= 18
y[contains.(lowercase.(meting),"_asb 65")] .= 18
y[contains.(lowercase.(meting),"_asb 66")] .= 18
y[contains.(lowercase.(meting),"_asb 67")] .= 18
y[contains.(lowercase.(meting),"_asb 68")] .= 18
y[contains.(lowercase.(meting),"_asb 69")] .= 19
y[contains.(lowercase.(meting),"_asb 70")] .= 19
y[contains.(lowercase.(meting),"_asb 71")] .= 19
y[contains.(lowercase.(meting),"_asb 72")] .= 20
y[contains.(lowercase.(meting),"_asb 73")] .= 20
y[contains.(lowercase.(meting),"_asb 74")] .= 20
y[contains.(lowercase.(meting),"_asb 75")] .= 19
y[contains.(lowercase.(meting),"_asb 76")] .= 19
y[contains.(lowercase.(meting),"_asb 77")] .= 19
y[contains.(lowercase.(meting),"_asb 78")] .= 21
y[contains.(lowercase.(meting),"_asb 79")] .= 21
y[contains.(lowercase.(meting),"_asb 80")] .= 21
y[contains.(lowercase.(meting),"_asb 81")] .= 8
y[contains.(lowercase.(meting),"_asb 82")] .= 8
y[contains.(lowercase.(meting),"_asb 83")] .= 8
y[contains.(lowercase.(meting),"_asb 84")] .= 8
y[contains.(lowercase.(meting),"_asb 85")] .= 8
y[contains.(lowercase.(meting),"_asb 86")] .= 8
y[contains.(lowercase.(meting),"_asb 87")] .= 8
y[contains.(lowercase.(meting),"_asb 88")] .= 8
y[contains.(lowercase.(meting),"_asb 89")] .= 7
y[contains.(lowercase.(meting),"_asb 90")] .= 7
y[contains.(lowercase.(meting),"_asb 91")] .= 7
y[contains.(lowercase.(meting),"_asb 92")] .= 7
y[contains.(lowercase.(meting),"_asb 93")] .= 7
y[contains.(lowercase.(meting),"_asb 94")] .= 7
y[contains.(lowercase.(meting),"_asb 95")] .= 7
y[contains.(lowercase.(meting),"_asb 96")] .= 22
y[contains.(lowercase.(meting),"_asb 97")] .= 2
y[contains.(lowercase.(meting),"_asb 98")] .= 23
class = fill("",length(y))

#SS1
class[y.==1] .= "MME Strobilanthes cusia (Nees) Kuntze"
class[y.==2] .= "MME Wrightia laevis Hook.f."
class[y.==3] .= "MME Indigofera suffruticosa Mill."
class[y.==4] .= "MME Indigofera tinctoria L."
class[y.==5] .= "LPE Strobilanthes cusia (Nees) Kuntze"
class[y.==6] .= "LPE Wrightia laevis Hook.f."
class[y.==7] .= "ADD Strobilanthes cusia (Nees) Kuntze"
class[y.==8] .= "M Isatis tinctoria L."
class[y.==9] .= "ADD Indigofera tinctoria L."
#SS2
class[y.==10] .= "HallTex Isatis tinctoria L."
class[y.==11] .= "ACIAR Strobilanthes cusia (Nees) Kuntze"
class[y.==12] .= "WrB Persicaria tinctoria Spach."
class[y.==13] .= "WrB Trade "
#SS3 Austria
class[y.==14] .= "MME. Persicaria tinctoria Spach."
class[y.==15] .= "MME. Indigofera tinctoria L."
class[y.==16] .= "MME. Indigofera suffruticosa Mill."
class[y.==17] .= "MME. Indigofera arrecta Hochst."
class[y.==18] .= "sLPE Persicaria tinctoria Spach."
class[y.==19] .= "sLPE Indigofera tinctoria L."
class[y.==20] .= "sLPE Indigofera arrecta Hochst."
class[y.==21] .= "sLPE Isatis tinctoria L."
class[y.==22] .= "ADD historic sample - unknown"
class[y.==23] .= "LPE Indigofera tinctoria L."

y2 = fill("",size(meting))
#SS1
y2[contains.(class, "MME ")] .= "diamond"
y2[contains.(class, "LPE")] .= "square"
y2[contains.(class, "ADD")] .= "star5"
y2[contains.(class, "M Isatis tinctoria L.")] .= "dtriangle"
#SS2
y2[contains.(class, "HallTex")] .= "utriangle"
y2[contains.(class, "ACIAR")] .= "pentagon"
y2[contains.(class, "WrB")] .= "hexagon"
y2[contains.(class, "WrB Trade")] .= "star6"
#SS3
y2[contains.(class, "MME.")] .= "circle"
y2[contains.(class, "sLPE")] .= "star4"
Symbol.(y2)

y3 = fill("",size(meting))
y3[contains.(class, "Strobilanthes cusia (Nees) Kuntze")] .= "medium purple2"
y3[contains.(class, "Wrightia laevis Hook.f.")] .= "mediumseagreen"
y3[contains.(class, "Indigofera suffruticosa Mill.")] .= "red"
y3[contains.(class, "Indigofera tinctoria L.")] .= "dodgerblue"
y3[contains.(class, "Isatis tinctoria L.")] .= "gold1"
y3[contains.(class, "Persicaria tinctoria Spach.")] .= "hot pink"
y3[contains.(class, "Indigofera arrecta Hochst.")] .= "palegreen"
y3[contains.(class, "HallTex Isatis tinctoria L.")] .= "gold1"
y3[contains.(class, "WrB Persicaria tinctoria Spach.")] .= "hot pink"
y3[contains.(class, "WrB Trade")] .= "grey0"
y3[contains.(class, "ADD historic sample - unknown")] .= "grey"
y3[contains.(class, "Isatis tinctoria L.")] .= "gold1"
Symbol.(y3)

## plotting                                              
PCA12 =
    scatter(scores[:,1], scores[:,2], group = class, markershape = Symbol.(y2), markercolor = Symbol.(y3), frame =:box, dpi = 120)
    scatter!(legend =:outerright, legendfont =8)
    #scatter!(legend=false)
    plot!([0], seriestype="vline", linecolor = :black, label="")
    plot!([0], seriestype="hline", linecolor = :black, label="")
    scatter!(xlabel= "PC 1 ($PC1"*"%)", ylabel="PC 2 ($PC2"*"%)", tickfont=font(12),guidefont=font(12))
    scatter!(title = "$Sample $Charge PCA 2D Plot - PC1 vs PC2", titlefont = 14)
    savefig("/Users/jianihu/Documents/PCA full negative/$Sample $Charge $plotv PCA_1vs2.pdf")

PCA13 = 
    scatter(scores[:,1], scores[:,3], group = class, markershape = Symbol.(y2), markercolor = Symbol.(y3), frame =:box, dpi = 120)
    scatter!(legend =:topleft, legendfont =7)
    #scatter!(legend = false)
    plot!([0], seriestype="vline", linecolor = :black, label="")
    plot!([0], seriestype="hline", linecolor = :black, label="")
    scatter!(xlabel= "PC 1 ($PC1"*"%)", ylabel="PC 3 ($PC3"*"%)", tickfont=font(12),guidefont=font(12))
    scatter!(title = "$Sample $Charge PCA 2D Plot - PC1 vs PC3", titlefont = 14)
    savefig("/Users/jianihu/Documents/PCA full negative/$Sample $Charge $plotv PCA_1vs3.pdf")

PCA14
    scatter(scores[:,1], scores[:,4], group = class, markershape = Symbol.(y2), markercolor = Symbol.(y3), frame =:box, dpi = 120)
    scatter!(legend = :topleft, legendfont = 6)
    #scatter!(legend = false)
    plot!([0], seriestype="vline", linecolor = :black, label="")
    plot!([0], seriestype="hline", linecolor = :black, label="")
    scatter!(xlabel= "PC 1 ($PC1"*"%)", ylabel="PC 4 ($PC4"*"%)", tickfont=font(12),guidefont=font(12))
    scatter!(title = "$Sample $Charge PCA 2D Plot - PC1 vs PC4", titlefont = 14)
    savefig("/Users/jianihu/Documents/PCA full negative/$Sample $Charge $plotv PCA_1vs4.pdf")

PCA23
    scatter(scores[:,2], scores[:,3], group = class, markershape = Symbol.(y2), markercolor = Symbol.(y3), frame =:box, dpi = 120)
    scatter!(legend = :topright, legendfont = 6)
    #scatter!(legend = false)
    plot!([0], seriestype="vline", linecolor = :black, label="")
    plot!([0], seriestype="hline", linecolor = :black, label="")
    scatter!(xlabel= "PC 2 ($PC2"*"%)", ylabel="PC 3 ($PC3"*"%)", tickfont=font(12),guidefont=font(12))
    scatter!(title = "$Sample $Charge PCA 2D Plot - PC2 vs PC3", titlefont = 14)
    savefig("/Users/jianihu/Documents/PCA full negative/$Sample $Charge $plotv PCA_2vs3.pdf")

PCA24
    scatter(scores[:,2], scores[:,4], group = class, markershape = Symbol.(y2), markercolor = Symbol.(y3), frame =:box, dpi = 120)
    scatter!(legend = :bottomleft, legendfont = 6)
    #scatter!(legend = false)
    plot!([0], seriestype="vline", linecolor = :black, label="")
    plot!([0], seriestype="hline", linecolor = :black, label="")
    scatter!(xlabel= "PC 2 ($PC2"*"%)", ylabel="PC 4 ($PC4"*"%)", tickfont=font(12),guidefont=font(12))
    scatter!(title = "$Sample $Charge PCA 2D Plot - PC2 vs PC4", titlefont = 14)
    savefig("/Users/jianihu/Documents/PCA full negative/$Sample $Charge $plotv PCA_2vs4.pdf")

PCA34
    scatter(scores[:,3], scores[:,4], group = class, markershape = Symbol.(y2), markercolor = Symbol.(y3), frame =:box, dpi = 120)
    scatter!(legend = :topright, legendfont = 5)    
    #scatter!(legend = false)
    plot!([0], seriestype="vline", linecolor = :black, label="")
    plot!([0], seriestype="hline", linecolor = :black, label="")
    scatter!(xlabel= "PC 3 (15.60%)", ylabel="PC 4 ($PC4"*"%)", tickfont=font(12),guidefont=font(12))
    scatter!(title = "$Sample $Charge PCA 2D Plot - PC3 vs PC4", titlefont = 14)
    savefig("/Users/jianihu/Documents/PCA full negative/$Sample $Charge $plotv PCA_3vs4.pdf")

## Extract loadings ## 
LoadPC1 = sortperm(abs.(loading[:, 1]))[end-9:end]
LoadPC1 = z[LoadPC1,:]
LPC1 = CSV.write("/Users/jianihu/Documents/PCA full negative/$T1 $Sample $Charge $plotv loading1.csv",LoadPC1)

LoadPC2 = sortperm(abs.(loading[:, 2]))[end-9:end]
LoadPC2 = z[LoadPC2,:]
LPC2 = CSV.write("/Users/jianihu/Documents/PCA full negative/$T1 $Sample $Charge $plotv loading2.csv",LoadPC2)

LoadPC3 = sortperm(abs.(loading[:, 3]))[end-9:end]
LoadPC3 = z[LoadPC3,:]
LPC3 = CSV.write("/Users/jianihu/Documents/PCA full negative/$T1 $Sample $Charge $plotv loading3.csv",LoadPC3)

LoadPC4 = sortperm(abs.(loading[:, 4]))[end-9:end]
LoadPC4 = z[LoadPC4,:]
LPC4 = CSV.write("/Users/jianihu/Documents/PCA full negative/$T1 $Sample $Charge $plotv loading4.csv",LoadPC4)
