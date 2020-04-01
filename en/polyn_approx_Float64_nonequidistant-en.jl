# tested for Julia v1.2
##### this version works with Float64 arithmetic
##### !!!! and for data with a (slightly) non-uniform time step (for the discrete argument)

# © Siberian Federal University, 2017
# Authors: Sergey P. Tsarev, Alexander S. Pustoshilov
# Paper reference: S.P.Tsarev, A.A.Kytmanov, Discrete orthogonal polynomials as a tool for detection of small anomalies of time series: a case study of GPS final orbits. arxiv.org/????? 2020

# This script robustly approximates a (probably vector-valued) time series with polynomials of large degrees
# (up to a few hundred) in the sense of the best RMS (best least squares) approximation
# any one-dimensional set of scalar or vector-valued data (which formally can be considered a function of a discrete argument)
# with a very large number of points (thousands)
# and find the approximation residue.
# You can approximate the entire dataset by polynomials on the entire data (or a part thereof) available in the file,
# or alternatively using a sliding window (possibly with overlapping).

# Requires auxiliary script service_functions-en.jl with some service functions

############################################################
# Script launch parameters:
############################################################

# julia.exe <script-name.jl> <file> <read_columns> <output_degrees> <windowsparams>
# where:

# <script-name.jl> - the name of this file
 
# <file> - name of the input file with the data being processed
# The file contains the data for approximation in columns - a plain text file with numeric data in any standard format: integer, floating point, exponential format.
# The file may contain empty and commented lines (the # character is used for commenting), such lines will be skipped in processing.
 
# <read_columns> - numbers of columns read from the file;
# numbers are indicated in square brackets,
# the first is the number of the column containing the discrete argument ("time"),
# if the line number should be used as a discrete argument,
# you have to specify 0 as the argument column (see examples);
# second and subsequent column numbers are columns that contain (vector)-function of the disctere argument.
# Examples of specifying column numbers:
# [1,2] - the first column of the file will be read as a discrete argument, the second - as a function of this argument;
# [0,1] - the line number will act as a discrete argument, the first column of the file as its function;
# [2,5:9] - the second column will act as a discrete argument, the  vector-function columns will be 5,6,7,8,9

# <output_degrees> - indicates degrees of approximating polynomials that will be used to output the results.
# Best least squares approximation is performed from degree 0 to the maximum degree specified in this parameter, but only the specified degrees will produce the output files.
# Examples of indicating the degrees:
# [1,2,3,4,5] - approximation by polynomials up to the 5th degree is performed, and the result is output after approximation by polynomials of degrees (1,2,3,4,5);
# [5] - approximation is carried out again by polynomials up to the 5th degree, but the result is output only for approximation by the best approximation polynomial of the 5th degree;
# [1:4, 5:5:20] - approximation by polynomials up to the 20th degree is performed, and the result is output for approximations by polynomials of degrees (1,2,3,4,5,10,15,20).
# The notations 1:4 and 5:5:20 have the same meaning as for loops in Julia.

# <windowsparams> - data processing windows parameters
# The script has the ability to process data in three ways:
# M:N - one window is given;
# M:W:N - windows with one point overlapping;
# M:W:S:N - with windows with a given overlap.
# Here:
# M - starting point of the first window,
# N - last point of the last window (as the last point
# you can specify the end of the file, the "end" keyword is used for this),
# W - window size
# S - windows overlap.
# !!! IMPORTANT NOTE - 1 !!!
# When forming windows, the last window is formed so that
# the last point of this window coincided with the last given point N,
# at the same time, the overlap with the penultimate window may be greater than specified in the parameter S.
# !!! IMPORTANT NOTE - 2 !!!
# If the specified window size is larger than the amount of data available, then the script exits logging the error in the log file.
# *************************************************** *********
# Examples of specifying window parameters:
# 1:end - one window is formed, the entire input file will be processed;
# 1:1000 - one window with a size of 1000 points is formed, the starting point;
# matches the first data point, the window end matches 1000th data point;
# 1:100:1000 - windows with the size of 100 points are formed at the following points (1, 100), (100, 199), (199, 298), ..., (901, 1000);
# 1:100:10:1000 - windows with a size of 100 points are formed with 10 points overlapping: (1, 100), (91, 190), (181, 280), ..., (901, 1000).

####################################################### ###########
# SCRIPT RESULTS
####################################################### ###########

# As a result of the script, APPROX directory is formed in the starting directory,
# in which one subdirectory is created with the name <out_dir_name>,
# then in this subdirectory individually for each approximated column from <read_columns>,
# and each degree from <output_degrees> and each processing window,
# a file containing 4 columns is created. Its columns are:
# <discrete arguments>,
# <initial approximated values>,
# <values ​​of the approximating polynomial>,
# <approximation residue = difference between the second and third columns>.
# The name of each file containd the number of the approximated column, the approximation window (line number of the file,
# not including comments) and the degree of approximating polynomial.
# All results in the files are printed with @sprintf ("% 1.15e", ...) - you can change the number of printed digits if necessary.
# For convenience of plotting the results using gnuplot we set the file extension: *.gplot_4
# *************************************************** *********
# Additionally, a file is created in the APPROX directory (appending new lines if the file exists) result_max_resid.txt,
# which contains the maximum(abs(residues)) at the maximal degree approximation for each window for each of the input columns,
# and the following data is written as one new line:
# - <out_dir_name> name of the output approximation directory;
# - <current_window> window being processed;
# - maximum(abs(residues)) for each approximated column in the current window with the maximal degree of approximation.

# *************************************************** *********
# To keep the errors, the errorlog_approx.log file is created in the script launch directory,
# which displays the following types of errors with the name of the experiment and the start time:
# "Incorrect data, 0 lines read" - occurs if the script could not
# read the source file or invalid column numbers were specified;
# "Invalid data, window length greater than data size" -
# informs that when setting windows parameters the window size was incorrectly selected
# and there is not enough data to cover at least 1 window.

####################################################### ###########
# Launch example:
####################################################### ###########
# julia <script.jl> F-testdata.txt [1:4] [0,1,3,5:5:20] 1:end
# 5:5:20 in <output_degrees> means a list of degrees with start 5, step 5 and end 20
# 1:end means processing the entire "time" interval, from the beginning of the file to the end as a single window.
####################################################### ###########

############################################################
############################################################

if length(ARGS) != 4
    println(length(ARGS))
    println("warning: wrong number of script parameters")
    exit()
end

using DelimitedFiles
using LinearAlgebra
using Printf
using Dates

# Loading the module with service functions

include("service_functions-en.jl")

# Beginning of the main approximation script

out_dir_name = "nonequ_$(ARGS[1])_col$(ARGS[2])_degs$(ARGS[3])_w[$(ARGS[4])]_Float64"
out_dir_name = replace(out_dir_name, ":"=>"-")

errorlog = "errorlog_approx.log"
max_resid_file = "APPROX/result_max_resid.txt"

filename = ARGS[1]
read_columns = make_numberlist(ARGS[2])

DATA = readbigfile(filename, read_columns, Float64)
length_DATA = length(DATA[:,1])

if length_DATA == 0
    add_log_messange(errorlog, "$(out_dir_name) Incorrect data, 0 lines read")
    exit()
end


output_degrees = make_numberlist(ARGS[3])
windows_list = make_windows(ARGS[4], length_DATA)

if length(windows_list) == 0
    add_log_messange(errorlog, "$(out_dir_name) Incorrect data, window size larger then the data size")
    exit()
end

mkpath("APPROX/$out_dir_name")

out_approx_max_resid = open(max_resid_file, "a")

cd("APPROX/$out_dir_name")

# maxdeg = output_degrees[end]
maxdeg = maximum(output_degrees)
sym_num_max_deg = length("$maxdeg")

# orthogonality tolerance: 
prec10 = 15
# number of binary digits 
prec = round(Int, prec10*log2(10))
orttol = 2.0^(-prec*90/100)
println("orttol=", orttol)

# Processing the sliding windows

for current_window in windows_list
    println("Processing window: $current_window")
    flush(stdout)
    # Preparing the table of discrete orthogonal polynomial values for the CURRENT window
    # pnts array contains the values of the discrete argument (the first specified column) of the one-dimensional discrete (vector) function (the remaining specified columns)
    pnts = DATA[current_window[1]:current_window[2],1]
    NN = length(pnts)
    p1 = pnts[1]
    pNN = pnts[end]
    # normalize the discrete argument to the interval [-1,1]:
    pnts11 = zeros(Float64, NN)
    for k=1:NN
        pnts11[k] = 2*(pnts[k] - p1) / (pNN - p1) - 1
    end
    pols = zeros(Float64, NN, maxdeg + 1)
    # pols[j, k] contains the value of orthogonal polynomial of degree (k-1) at the jth point of the discrete argument (specified by the first column)
    # At first we fill the array pols with the values of Legendre polynomials (they are NOT orthogonal on the discrete lattice!), then we use a special variant of Gram-Schmidt orthogonalization (and the we normalize the polynomials)
    # first polynomial is constant:
    pols[:,1] = ones(Float64, NN)
    # the second is a linear function of the discrete argument (specified by the first column)
    pols[:,2] = pnts11
    for d = 3:(maxdeg + 1)
        n = d - 2
        # pols[:,d] =   (((2 * n + 1) .* pnts11 .* pols[:,d - 1]) - (n .* pols[:,d - 2])) ./ (n + 1);
        pols[:,d] =   (((2 * n + 1) * (pnts11 .* pols[:,d - 1])) - (n * pols[:,d - 2])) / (n + 1);
        println("degree = ", d - 1)
        flush(stdout)
    end
    println("Legendre polynomials done")
    norm0 = sqrt(sum(abs2, pols[:,1]))
    pols[:,1] = pols[:,1] / norm0
    println("Orthogonalization started")
    flush(stdout)
    max_resid_in_col = zeros(Float64, length(read_columns)-1)
    # The obtained Legendre polynomials now need to be orthogonalized and normalized.
    # pols[:,2] SHOULD be orthogonalized for non-equidistant points  !!!!
    for n = 2:(maxdeg + 1)
        u = pols[:,n]
        repN = 0
        nev = 1000000.0
        while ( nev > 2 * NN * orttol )
            repN = repN + 1
            for i = 1:(n - 1)
                cf = sum(u .* pols[:,i])
                u = u - cf * pols[:,i]  # Wikipedia version
                # - this version of Gram-Schmidt algorithm turned out to be robust even for Float64 computations!!!!
            end
            nev = 0.0
            norm0 = sqrt(sum(abs2, u))
            for i = 1:(n - 1)
                # nev = max(nev, abs(sum(u .* pols[:,i])))
                nev = max(nev, abs(sum(u .* pols[:,i]))/norm0)
            end
        end
        norm0 = sqrt(sum(abs2, u))
        pols[:,n] = u / norm0
        println("degree = ", n - 1, "   max nev after orthogonalization = ", @sprintf("%1.15e",nev),
        "   repN = ", repN, "   norm was: ", @sprintf("%1.15e",norm0))
        flush(stdout)
    end
    # End of orthogonalization - in the array pols[:,:] we have the set of values of the normalized discrete orthogonal Hahn polynomials for the CURRENT window
    println("Orthogonalization done")
    flush(stdout)
    for ind_cn = 2:length(read_columns)
        data_col = DATA[current_window[1]:current_window[2],ind_cn]
        data_arg = DATA[current_window[1]:current_window[2], 1]
        data_cn = copy(data_col)
        pol_approx = zeros(Float64, NN)
        # approximation of the input data (individually for each column) by polynomials - since they are orthonormal, a simple scalar product is enough!
        for dn in 0:maxdeg  
            p = pols[:,dn + 1]
            cf = sum(data_cn .* p)
            data_cn = data_cn .- (cf * p)
            pol_approx = pol_approx .+ (cf * p)
            print("col: ", read_columns[ind_cn], "   win: ", current_window, "   deg: ", dn)
            flush(stdout)
            if (dn in output_degrees) 
                output_filename = "col$(read_columns[ind_cn])_win$(current_window)_deg$(print_with_zero(dn, sym_num_max_deg)).gplot_4"
                output_filename = replace(output_filename, " "=>"_")
                output_filename = replace(output_filename, ","=>"")
                out_file = open(output_filename, "w")
                for k = 1:length(pnts)
                    @printf(out_file, "%1.15e\t%1.15e\t%1.15e\t%1.15e\n", data_arg[k], data_col[k], pol_approx[k], data_cn[k])
                end
                println(" written!")
                close(out_file)
            else
                println()
                flush(stdout)
            end            
        end # dn in 0:maxdeg
        max_resid_in_col[ind_cn-1] = maximum(abs.(data_cn))
    end # ind_cn
    println("$out_dir_name $(current_window) $(join(map(x->@sprintf("%1.15e",x),max_resid_in_col), "\t"))")   
    println(out_approx_max_resid, "$out_dir_name $(current_window) $(join(map(x->@sprintf("%1.15e",x),max_resid_in_col), "\t"))")   
    flush(out_approx_max_resid)            
end # current_window

close(out_approx_max_resid)
exit()
