## PLS DA
using ScikitLearn, Plots, CSV, DataFrames, Statistics
@sk_import cross_decomposition: PLSRegression

## Data import ## 
#neg 
data = CSV.read("/Users/jianihu/Documents/New featurelists/FeatureList_Aligned_SS3_neg_full_PLSDA.csv", DataFrame)

#pos
data = CSV.read("/Users/jianihu/Documents/New featurelists/PLSDA csv/FeatureList_Aligned_SS3_pos_full_PLSDA.csv", DataFrame)

b = Matrix(data)'

function standard(b)
    m = mean(b,dims=1)                          # mean berekenen
    stdev = std(b,dims=1)                       # stdev berekenen
    bs = (b .- m) ./ stdev                      # van alles de mean aftrekken en delen door de stdev => meancentering & scaling
    return bs
end

s_data = standard(b)

## classification & plot labeling ## 
meting = names(data)
y = zeros(size(meting))
y[contains.(lowercase.(meting),"_asb 24")] .= 1
y[contains.(lowercase.(meting),"_asb 25")] .= 1
y[contains.(lowercase.(meting),"_asb 26")] .= 1
y[contains.(lowercase.(meting),"_asb 27")] .= 1
y[contains.(lowercase.(meting),"_asb 28")] .= 1
y[contains.(lowercase.(meting),"_asb 29")] .= 1
y[contains.(lowercase.(meting),"_asb 30")] .= 1
y[contains.(lowercase.(meting),"_asb 31")] .= 1
y[contains.(lowercase.(meting),"_asb 32")] .= 1

y[contains.(lowercase.(meting),"_asb 33")] .= 2
y[contains.(lowercase.(meting),"_asb 34")] .= 2
y[contains.(lowercase.(meting),"_asb 35")] .= 2

y[contains.(lowercase.(meting),"_asb 36")] .= 3
y[contains.(lowercase.(meting),"_asb 37")] .= 3
y[contains.(lowercase.(meting),"_asb 38")] .= 3

y[contains.(lowercase.(meting),"_asb 39")] .= 2
y[contains.(lowercase.(meting),"_asb 40")] .= 2
y[contains.(lowercase.(meting),"_asb 41")] .= 2
y[contains.(lowercase.(meting),"_asb 42")] .= 2
y[contains.(lowercase.(meting),"_asb 43")] .= 2
y[contains.(lowercase.(meting),"_asb 44")] .= 2

y[contains.(lowercase.(meting),"_asb 45")] .= 3
y[contains.(lowercase.(meting),"_asb 46")] .= 3
y[contains.(lowercase.(meting),"_asb 47")] .= 3
y[contains.(lowercase.(meting),"_asb 48")] .= 3
y[contains.(lowercase.(meting),"_asb 49")] .= 3
y[contains.(lowercase.(meting),"_asb 50")] .= 3

y[contains.(lowercase.(meting),"_asb 51")] .= 4
y[contains.(lowercase.(meting),"_asb 52")] .= 4
y[contains.(lowercase.(meting),"_asb 53")] .= 4
y[contains.(lowercase.(meting),"_asb 54")] .= 4
y[contains.(lowercase.(meting),"_asb 55")] .= 4
y[contains.(lowercase.(meting),"_asb 56")] .= 4

y[contains.(lowercase.(meting),"_asb 57")] .= 2
y[contains.(lowercase.(meting),"_asb 58")] .= 2
y[contains.(lowercase.(meting),"_asb 59")] .= 2

y[contains.(lowercase.(meting),"_asb 60")] .= 5
y[contains.(lowercase.(meting),"_asb 61")] .= 5
y[contains.(lowercase.(meting),"_asb 62")] .= 5
y[contains.(lowercase.(meting),"_asb 63")] .= 5
y[contains.(lowercase.(meting),"_asb 64")] .= 5
y[contains.(lowercase.(meting),"_asb 65")] .= 5
y[contains.(lowercase.(meting),"_asb 66")] .= 5
y[contains.(lowercase.(meting),"_asb 67")] .= 5
y[contains.(lowercase.(meting),"_asb 68")] .= 5

y[contains.(lowercase.(meting),"_asb 69")] .= 6
y[contains.(lowercase.(meting),"_asb 70")] .= 6
y[contains.(lowercase.(meting),"_asb 71")] .= 6

y[contains.(lowercase.(meting),"_asb 72")] .= 7
y[contains.(lowercase.(meting),"_asb 73")] .= 7
y[contains.(lowercase.(meting),"_asb 74")] .= 7

y[contains.(lowercase.(meting),"_asb 75")] .= 6
y[contains.(lowercase.(meting),"_asb 76")] .= 6
y[contains.(lowercase.(meting),"_asb 77")] .= 6

y[contains.(lowercase.(meting),"_asb 78")] .= 8
y[contains.(lowercase.(meting),"_asb 79")] .= 8
y[contains.(lowercase.(meting),"_asb 80")] .= 8

class = fill("",length(y))

class[y.==1] .= "MME Persicaria tinctoria"
class[y.==2] .= "MME Indigofera tinctoria"
class[y.==3] .= "MME Indigofera suffruticosa"
class[y.==4] .= "MME Indigofera arrecta"

class[y.==5] .= "sLPE Persicaria tinctoria"
class[y.==6] .= "sLPE Indigofera tinctoria"
class[y.==7] .= "sLPE Indigofera arrecta"
class[y.==8] .= "sLPE Isatis tinctoria"

y2 = fill("",size(meting))
y2[contains.(class, "MME")] .= "circle"
y2[contains.(class, "sLPE")] .= "start"
y2[contains.(class, "ADD")] .= "diamond"
Symbol.(y2)

y3 = fill("",size(meting))
y3[contains.(class, "Woad")] .= "blue"
y3[contains.(class, "Strobilanthes cusia")] .= "yellow"
y3[contains.(class, "Persicaria tinctoria Spach.")] .= "red"
y3[contains.(class, "trade")] .= "green"
Symbol.(y3)

## create binaire matrix ## 
y_pls = zeros(size(data,2),8)
for i = 1:8
    y_pls[y.==i,i] .= 1
end

PLSR = PLSRegression(n_components=8)

b[b .> 0].= 1

## PLS-DA code ##
Y = y_pls

PLSR.fit(b, Y)

PLSR.predict(b)

y_s = PLSR.y_scores_

## kwantitatief -
scatter(y_s[:,1],y_s[:,2],group=class, frame =:box, dpi = 120)
scatter!(xlabel= "y scores 1", ylabel="y scores 2", legend = :topleft, legendfont = 6, frame = :box, dpi =120) #legend =:bottomleft, legendfont = 7
scatter!(title = "PLS-DA 2D plot - SS3 Quantitative negative")
savefig("/Users/jianihu/Documents/SS2_negatief_PLSDA_2D_kwanti_1vs2.pdf")

scatter(y_s[:,1],y_s[:,2],y_s[:,3], group=class, frame =:box, dpi = 120,markershape = Symbol.(y2), markercolor = Symbol.(y3))
scatter!(xlabel= "y scores 1", ylabel="y scores 2", zlabel = "y scores 3", legend =:bottomright, legendfont = 7, frame = :box, dpi =120) #legend =:bottomleft, legendfont = 7
scatter!(title = "PLS-DA 3D plot - SS2 Quantitative negative")
savefig("/Users/jianihu/Documents/SS2_negatief_PLSDA_3D_kwanti_1vs2vs3.pdf")

## kwantitatief +
scatter(y_s[:,1],y_s[:,2],group=class, frame =:box, dpi = 120, markershape = Symbol.(y2), markercolor = Symbol.(y3))
scatter!(xlabel= "y scores 1", ylabel="y scores 2", legend = true, legendfont = 8, frame = :box, dpi =120) #legend =:bottomleft, legendfont = 7
scatter!(title = "PLS-DA 2D plot - SS2 Quantitative positive")
savefig("/Users/jianihu/Documents/SS2_positief_PLSDA_2D_kwanti_1vs2.pdf")

scatter(y_s[:,1],y_s[:,2],y_s[:,3], group=class, frame =:box, dpi = 120,markershape = Symbol.(y2), markercolor = Symbol.(y3))
scatter!(xlabel= "y scores 1", ylabel="y scores 2", zlabel = "y scores 3", legend =:topright, legendfont = 8, frame = :box, dpi =120) #legend =:bottomleft, legendfont = 7
scatter!(title = "PLS-DA 3D plot - SS2 Quantitative positive")
savefig("/Users/jianihu/Documents/SS2_positief_PLSDA_3D_kwanti_1vs2vs3.pdf")

## kwalitatief -
scatter(y_s[:,1],y_s[:,2],group=class, frame =:box, dpi = 120)
scatter!(xlabel= "y scores 1", ylabel="y scores 2", legend = false, frame = :box, dpi =120) #legend =:bottomleft, legendfont = 7
scatter!(title = "PLS-DA 2D plot - SS1 Qualitative negative")
savefig("/Users/jianihu/Documents/SS1_negatief_PLSDA_2D_kwanli_1vs2.pdf")

scatter(y_s[:,1],y_s[:,2],y_s[:,3], group=class, frame =:box, dpi = 120)
scatter!(xlabel= "y scores 1", ylabel="y scores 2", zlabel = "y scores 3", legend = false, frame = :box, dpi =120) #legend =:bottomleft, legendfont = 7
scatter!(title = "PLS-DA 3D plot - SS1 Qualitative negative")
savefig("/Users/jianihu/Documents/SS1_negatiefPLSDA_3D_kwali_1vs2vs3.pdf")
## kwalitatief +
scatter(y_s[:,1],y_s[:,2],group=class, frame =:box, dpi = 120)
scatter!(xlabel= "y scores 1", ylabel="y scores 2", legend = false, frame = :box, dpi =120) #legend =:bottomleft, legendfont = 7
scatter!(title = "PLS-DA 2D plot - SS1 Qualitative positive")
savefig("/Users/jianihu/Documents/SS1_positief_PLSDA_2D_kwanli_1vs2.pdf")

scatter(y_s[:,1],y_s[:,2],y_s[:,3], group=class, frame =:box, dpi = 120)
scatter!(xlabel= "y scores 1", ylabel="y scores 2", zlabel = "y scores 3", legend = false, frame = :box, dpi =120) #legend =:bottomleft, legendfont = 7
scatter!(title = "PLS-DA 3D plot - SS1 Qualitative positive")
savefig("/Users/jianihu/Documents/SS1_positief_PLSDA_3D_kwali_1vs2vs3.pdf")

## PLS-DA
PLSR.fit_transform(b, Y)

xw = PLSR.x_weights_ #The left singular vectors of the cross-covariance matrices of each iteration
yw = PLSR.y_weights_ #The right singular vectors of the cross-covariance matrices of each iteration

xs = PLSR.x_scores_
ys = PLSR.y_scores_

xl = PLSR.x_loadings_
yl = PLSR.y_loadings_

## Extract loadings
LPLSDA1 = sortperm(abs.(PLSR.x_loadings_[:, 1]))[end-9:end]
LPLSDA1 = z[LPLSDA1,:]
PLSDA1 = CSV.write("/Users/jianihu/Documents/New featurelists/PLSDA SS2 neg loading1 +.csv",LPLSDA1)

LPLSDA2 = sortperm(abs.(PLSR.x_loadings_[:, 2]))[end-9:end]
LPLSDA2 = z[LPLSDA2,:]
PLSDA2 = CSV.write("/Users/jianihu/Documents/New featurelists/PLSDA SS2 neg loading2 +.csv",LPLSDA2)

LPLSDA3 = sortperm(abs.(PLSR.x_loadings_[:, 3]))[end-9:end]
LPLSDA3 = z[LPLSDA3,:]
PLSDA3= CSV.write("/Users/jianihu/Documents/New featurelists/PLSDA SS2 neg loading3 +.csv",LPLSDA3)

LPLSDA4 = sortperm(abs.(PLSR.x_loadings_[:, 4]))[end-9:end]
LPLSDA4 = z[LPLSDA4,:]
PLSDA4 = CSV.write("/Users/jianihu/Documents/New featurelists/PLSDA SS2 neg loading4+.csv",LPLSDA4)
