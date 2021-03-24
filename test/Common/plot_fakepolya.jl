using Plots, Printf
import DelimitedFiles

pyplot() # plotting backend

include("plot_defs.jl")

# plot the data for the simple periodic example

dirname = @__DIR__
ph = plot()
for i = 4:2:16
    fname = @sprintf "%s%s%03d%s" dirname "/output/fakepolya" i ".txt" 
    aa = DelimitedFiles.readdlm(fname)
    plotlab = @sprintf "Order %d" i
    plot!(aa[2:end,3] ./ aa[1,4], aa[2:end,6], label=plotlab)
end
plot!([1.25,2],[5,5],label="5",style=:dash,linewidth=2)

xaxis!("number of pts/number of polys")
yaxis!("normalized condition number",:log10)
title!("stability vs efficiency of thinned nodes")

#  save figure
mkpath(dirname * "/figures")
fname = dirname * "/figures/fakepolyastab.pdf"
savefig(ph,fname)
