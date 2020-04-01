# This file contains utility functions necessary for:
# 1. handling input arguments to a script call
# 2. read large files in a given format with possible omission of commented and unnnecessary initial lines
# 3. output a number of type Int with the addition of the required number of zeros at the beginning
# 4. event logging function

# Functions for processing an incoming script argument specifying the degrees of approximation for output
# This function processes the parameter enclosed in square brackets
# and based on the specified data forms an array of elements
# Example:
# if the string "[1,2,3,4,5: 5: 20,30,40]" is input,
# Array{Int}[1,2,3,4,5,10,15,20,30,40] is formed as the output
function make_numberlist(arg_string)
    
    arg_split = split(strip(arg_string, ['[',']']), ",")

    arg_deg = Any[]
    for arg in arg_split
        s_arg = split(arg, ":")
        if length(s_arg) >1
            if length(s_arg) == 2
                st = parse(Int, s_arg[1])
                en = parse(Int, s_arg[2])
                push!(arg_deg, st:en)
            end
            if length(s_arg) == 3
                st = parse(Int, s_arg[1])
                sp = parse(Int, s_arg[2])
                en = parse(Int, s_arg[3])
                push!(arg_deg, st:sp:en)
            end
        else
            push!(arg_deg, parse(Int, arg))
        end
    end
    degs = Int[]
    for els in arg_deg
        if typeof(els) <: Integer
            push!(degs, els)
        else
            for el_in_collect in els
                push!(degs, el_in_collect)
            end
        end
    end
    return degs
end

# This function is designed to form sliding windows
# The input of the function is a string of arguments with windows parameters and the maximum size of the data onto which the windows are formed
# The argument string can be specified in three ways:
# M:N - in one window;
# M:W:N - windows overlapping in one point;
# M:W:S:N - with windows with a given overlap.
# Here:
# M - starting point of the first window
# N - last point of the last window (as the last point
# you can specify the end of the file, the "end" keyword is used for this)
# W - window size
# S - windows overlap
# !!! IMPORTANT NOTE - 1 !!!
# When forming windows, the last window is formed so that
# the last point of the window coincided with the last given point N,
# at the same time, the overlap with the penultimate window may be greater than specified S
# !!! IMPORTANT NOTE - 2 !!!
# If the specified window size is larger than the amount of data available, then the script exits logging the error in the log file.
# ***************************************************
# Examples
# "1:10", 20 => [(1,10)]
# "1:10:19", 18 => [(1, 10), (9, 18)]
# "1:10:19", 20 => [(1, 10), (10, 19)]
# "1:10:20", 20 => [(1, 10), (10, 19), (11, 20)]
# "1: 10: 3: 20", 20 => [(1, 10), (8, 17), (11, 20)]

function make_windows(arg_string, data_length)
    windows_array = Any[]
    arg = split(arg_string, ':')
    if length(arg) == 2 # Option when 1 window is used
        start_index = parse(Int, arg[1])
        if arg[2] == "end"
            end_index = data_length
        else
            end_index =  parse(Int, arg[2])
            if end_index>data_length
                return windows_array
            end
        end
        windows_array = [(start_index, end_index)]
    elseif length(arg) == 3 # windows with 1 point overlap
        start_index = parse(Int, arg[1])
        step_size = parse(Int, arg[2])
        if arg[3] == "end"
            end_index = data_length
        else
            end_index =  parse(Int, arg[3])
            if end_index>data_length
                return windows_array
            end            
        end
        for start_window = start_index:step_size-1:end_index - step_size        
            push!(windows_array, (start_window, start_window + step_size-1))
        end
        if length(windows_array) ==0
            return windows_array
        end
        if windows_array[end][2] < end_index
            push!(windows_array, (end_index - step_size+1, end_index))
        end
    elseif length(arg) == 4 # windows with overlapping
        start_index = parse(Int, arg[1])
        step_size = parse(Int, arg[2])
        overlap = parse(Int, arg[3])
        if arg[4] == "end"
            end_index = data_length
        else
            end_index =  parse(Int, arg[4])
            if end_index>data_length
                println("return")
                return windows_array
            end            
        end
        println("read param: start$start_index - step$step_size - overlap$overlap - end$end_index")
        for start_window = start_index:(step_size - overlap):(end_index - step_size)
            push!(windows_array, (start_window, start_window + step_size-1))
        end
        if length(windows_array) ==0
            return windows_array
        end
        if windows_array[end][2] < end_index
            push!(windows_array, (end_index - step_size + 1, end_index))
        end
    else
        windows_array = [(1, data_length)]
    end
    return windows_array
end

# function for data reading
function readbigfile(filename::AbstractString, # input file name
    Columns::Array{Int, 1}, # numbers of columns to be read. If you specify 0 as one of the numbers to be read, then a column with line numbers will be generated (lines omitted during reading are not considered)
    T::Type # data type into which columns will be converted
    ;
    skip_first_row = 0 # the number of lines at the beginning of the file to skip
    )::Array{T, 2} # function returns a two-dimensional array

    # Open the file for reading
    in_file = readlines(filename)
    skip_row = 0
    for line in in_file
        if skip_first_row > skip_row
            skip_row=skip_row+1
            continue
        end
        # skip lines starting with '#' (IN THE FIRST POSITION!)
        if line[1] == '#'
            skip_row = skip_row+1
            continue
        end
        # if length(split(line)) < length(Columns)
        #     skip_row = skip_row+1
        #     continue
        # end
    end
    println("skipped rows: $skip_row")
    # Create a zero two-dimensional array of type T
    out_array = zeros(T, (length(in_file) - skip_row, length(Columns)))
    # We process the file line by line
    num_out_array_row = 1
    for num_str = 1+skip_first_row:length(in_file)
        if in_file[num_str][1] == '#' 
            continue
        end
# skip rows with insufficient number of columns
        # if length(split(in_file[num_str])) < length(Columns)
        #     continue
        # end
        value_str = split(in_file[num_str])
        for num_col = 1:length(Columns)
            if Columns[num_col] == 0
                out_array[num_out_array_row,num_col] = num_out_array_row
            else    
                out_array[num_out_array_row,num_col] = parse(T, value_str[Columns[num_col]])
            end
        end
        num_out_array_row = num_out_array_row + 1
    end
    # Return the array of data
    return out_array
end

# This function takes two numbers as input:
# 1. number to be printed
# 2. minimum number of characters in the number to be displayed (as a string). If the printed form of the number (first argument) is
# shorter than specified, zeros are padded before the number so that the length of the printed number with zeros is equal to the specified minimum length
# If the number is longer, just the number itself is output as a string
function print_with_zero(num::Integer, min_length_num::Integer)
    str_num = "$num"    
    zero_in_begin = "0"^maximum([0,(min_length_num-length(str_num))])
    return "$(zero_in_begin)$(str_num)"
end


# The function appends the message indicated as the second argument to a file with the name specified as the first argument of the function,
# the current time stamp is indicated as well
function add_log_messange(log_file, message)
    out_file = open(log_file, "a")
    println(out_file, "[$(now())] $(message)")
    flush(out_file)
    close(out_file)
end

# Function outputs (by default to stdout) two-dimensional and one-dimensional arrays with the given separator between elements
# The input array is input and separators between the elements of the string can be specified (col_delims)
# as well as separators between the lines (row_delims)
function println_array(InputArray::Array;col_delims = "\t", row_delims="\n", out_stream = stdout)
    dims = ndims(InputArray)
    if dims == 1
        println(out_stream,join(map(x->"$(x)",InputArray),col_delims))
    elseif dims == 2
        num_rows = size(InputArray)[1]
        for ind_row = 1:num_rows
            print(out_stream,join(map(x->"$(x)",InputArray[ind_row, :]),col_delims))            
            print(out_stream,row_delims)
        end            
        println(out_stream)
    else
        println("warning: Dims > 2")
    end
end
