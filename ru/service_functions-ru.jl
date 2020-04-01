# В данном файле собраны служебные функции, необходимые для:
# 1. обработки входных аргументов вызова скрипта;
# 2. чтения больших файлов в заданном формате с возможным пропуском закомментированных или ненужных строк  заголовка;
# 3. вывод числа типа Int с добавлением необходимого числа нулей в начале;
# 4. функция логирования событий.

# Функции обработки входящего аргумента скрипта, задающего набор степеней полиномов для вывода результатов аппроксимации.
# Данная функция обрабатывает параметр, заключенный в квадратные скобки 
# и на основании указанных данных формирует массив элементов 
# Пример: 
# на вход подается строка "[1,2,3,4,5:5:20,30,40]"
# на выходе формируется массив Array{Int}[1,2,3,4,5,10,15,20,30,40]

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

# Данная функция предназначена для формирования набора окон для обработки
# На вход функции подается строка аргументов с параметрами окна и максимальный размер данных, на которые формируются окна
# Строка аргументов может быть задана тремя способами:
# M:N - одним окном; 
# M:W:N - окнами с перекрытием в одну точку;
# M:W:S:N - окнами с заданным перекрытием.
# При этом всюду выше:
# M - начальная точка первого окна,
# N - последняя точка последнего окна (в качестве последней точки 
# можно указать конец файла, для этого используется ключевое слово end),
# W - размер окна,
# S - перекрытие соседних окон.
# !!!ВАЖНОЕ ЗАМЕЧАНИЕ - 1!!!
# При формировании окон последнее окно формируется так, чтобы 
# последняя точка окна совпадала с последней заданной точкой N,
# при этом перекрытие с предпоследним окном может быть больше, чем указано в исходных параметрах
# !!!ВАЖНОЕ ЗАМЕЧАНИЕ - 2!!!
# Если размер окна задан больше, чем количество доступных данных, то скрипт выдает ошибку в лог-файл и завершает работу
# **********************************************************# Примеры
# "1:10", 20 => [(1,10)]
# "1:10:19", 18 => [(1, 10), (9, 18)]
# "1:10:19", 20 => [(1, 10), (10, 19)]
# "1:10:20", 20 => [(1, 10), (10, 19), (11, 20)]
# "1:10:3:20", 20 => [(1, 10), (8, 17), (11, 20)]

function make_windows(arg_string, data_length)
    windows_array = Any[]
    arg = split(arg_string, ':')
    if length(arg) == 2 # Вариант, когда используется 1 окно
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
    elseif length(arg) == 3 # окна c перекрытием в 1 точку
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
    elseif length(arg) == 4 # окна с перекрытием
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

# функция для чтения данных
function readbigfile(filename::AbstractString,           # имя входного файла
    Columns::Array{Int,1},     # номера стобцов, которые необходимо считать. Если указать 0 как один из номеров для считывания, то будет сформирован столбец с номерами строк (пропущенные при чтении строки не считаются)
    T::Type                    # тип данных, в который будут сконвертированы считанные столбцы
    ;
    skip_first_row = 0        # число строк в начале файла, которые необходимо пропустить
    )::Array{T,2}             # функция возвращает двумерный массив

    # Открываем файл для чтения
    in_file = readlines(filename)
    skip_row = 0
    for line in in_file
        if skip_first_row > skip_row
            skip_row=skip_row+1
            continue
        end
        # пропускаем строки, начинающиеся (В ПЕРВОЙ ПОЗИЦИИ!) c #
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
    # Создаем нулевой двумерный массив типа T
    out_array = zeros(T, (length(in_file) - skip_row, length(Columns)))
    # Построчно обрабатываем файл
    num_out_array_row = 1
    for num_str = 1+skip_first_row:length(in_file)
        if in_file[num_str][1] == '#' 
            continue
        end
#  пропускаем строки с недостаточным числом столбцов
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
    # Возвращаем считанный массив
    return out_array
end

# Данная функция принимает на вход два числа: 
# 1. целое положительное число, которое необходимо напечатать
# 2. минимальное число символов в выводимом числе. 
# Если символов в числе меньше указанного, перед числом добавляются нули так, чтобы число с нулями стало необходимой минимальной длины.
# Если число длиннее, выводится просто указанное число
function print_with_zero(num::Integer, min_length_num::Integer)
    str_num = "$num"    
    zero_in_begin = "0"^maximum([0,(min_length_num-length(str_num))])
    return "$(zero_in_begin)$(str_num)"
end


# Функция дописывает сообщение, указанное вторым аргументом функции в файл с именем указанным как первый аргумент функции,
# при этом указывается время записи
function add_log_messange(log_name, message)
    out_file = open(log_name, "a")
    println(out_file, "[$(now())] $(message)")
    flush(out_file)
    close(out_file)
end

# Функция для вывода двумерных и одномерных массивов с указанием разделителя между элементами
# На вход подается входной массив и могут быть указаны разделители между элементами строки (col_delims) 
# и между строками (row_delims), 
# если не указывать - по умолчанию между элементами табуляция и между строками перенос строки
# также может быть указан поток вывода, по умолчанию STDOUT
function println_array(InputArray::Array;col_delims = "\t", row_delims="\n", out_stream = STDOUT)
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
