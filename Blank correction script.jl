##### Blank correction script #####

## needed packages
using DataFrames, XLSX, Statistics, CSV, LinearAlgebra, Plots

## import FeatureList
ssc_neg_Intensities = DataFrame(XLSX.readtable("/Users/jianihu/Documents/FeatureList_Aligned_Full_Neg.xlsx","Intensities")...)      # Data import
ssc_neg_Intensities = ssc_neg_Intensities[ 3 .< ssc_neg_Intensities[!,"AverageRT"] .<30 ,:]                                         # Selecting chromatogram data between 3.0 - 30 min

## file in sample information - Change per imput
Sample = "SS All"       # Sample set 
Charge = "Negative"     # Negative or Positive mode 
plotv = "V1"            # version of plots 
T1 = "PCA"              # PCA or PLSDA

## Matrix clean up
# Blanco signal correction
function blanco_clean_up(data_mat)                      # data_mat = ss1_pos_Intensities of een andere
    samples = zeros(Bool,size(data_mat,2))
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
    
    samples[contains.(names(data_mat),"SS2").& (.! blanco) .& (.! standaard)] .= 1 # First part selecting the data and second part erasing/converting the other data to 0 
    
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

    return data_mat, samples                            # return data_mat and samples

end

z,samples = blanco_clean_up(ssc_neg_Intensities)        # z,samples zorgt ervoor dat de  
z = z[vec(any(Matrix(z)[:,samples].!=0,dims=2)),:]      # alle rijen, die gelijk zijn aan -0 weg

## Save import data after blank correction 
SSC_neg = CSV.write("/Users/jianihu/Documents/PCA full negative/FeatureList_Aligned_$Charge _full_PLSDA $Sample $plotv.csv", z)

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
