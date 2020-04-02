# tested for Julia v1.2
##### данная версия работает с арифметикой Float64
##### !!!!!! и только с данными с РАВНОМЕРНЫМ шагом "меток времени" (дискретного аргумента)

# © Сибирский федеральный университет, 2017
# Авторы: С.П.Царев, А.С.Пустошилов
# Теоретические детали и пример применения: S.P.Tsarev, A.A.Kytmanov, Discrete orthogonal polynomials as a tool for detection of small anomalies of time series: a case study of GPS final orbits: https://arxiv.org/abs/2004.00414

# Данный скрипт позволяет аппроксимировать полиномами больших степеней (степенью до нескольких сотен) 
# одномерный набор скалярных и векторнозначных данных (который формально можно считать функцией дискретного аргумента, или временным рядом)
# в смысле наилучшего среднеквадратичного приближения с очень большим количеством точек (тысячи) 
# и находить невязку аппроксимации.
# При этом можно приближать как весь отрезок данных полиномами на всем имеющемся в файле отрезке данных (или его части),
# так и проходя по нему скользящим окном (возможно, с перекрытием)

# Требует наличия вспомогательного скрипта service_functions-ru.jl со служебными функциями

############################################################
# Параметры запуска скрипта:
############################################################

# julia.exe <script-name.jl> <file> <read_columns> <output_degrees> <windowsparams>
# где:

# <script-name.jl> - имя данного файла

# <file> - имя обрабатываемого файла с данными.
# В качестве обрабатываемого файла должен выступать текстовый файл, содержащий данные для аппроксимации в отдельных колонках, в любом из числовых форматов (целые, с плавающей точкой, либо в экспоненциальном формате).
# Файл может содержать пустые и закомментированные строки (для комментирования используется символ #), такие строки в обработке будут пропущены.

# <read_columns> - номера считываемых из файла колонок 
# номера указываются в квадратных скобках
# первым указывается номер колонки, содержащей дискретный аргумент ("метка времени").
# Если в качестве дискретного аргумента должен выступать номер строки,
# то необходимо указать 0 как столбец аргумента (см. примеры).
# Вторым и последующими указываются номера колонок, содержащие (вектор)-функцию этого аргумента.
# Примеры указания номеров колонок:
# [1,2] - первый столбец файла будет считан как дискретный аргумент, второй - как скалярная функция от этого аргумента;
# [0,1] - в качестве дискретного аргумента будет выступать номер строки, в качестве функции - первый столбец файла;
# [2,5:9] - в качестве дискретного аргумента будет выступать второй столбец, в качестве вектор-функции - столбцы 5,6,7,8,9 

# <output_degrees> - какие степени полиномов будут использованы для вывода результатов.
# Фактически аппроксимация выполняется от 0 до максимальной степени, указанной в данном параметре, но выводится будут результаты только для указанных степеней.
# Примеры указания выводимых степеней:
# [1,2,4,5] - выполняется аппроксимация полиномами до 5-ой степени, и выводится результат после аппроксимации полиномами степени 1,2,4,5;
# [5] - выполняется аппроксимация полиномами до 5-ой степени, и выводится результат аппроксимации полиномами 5-ой степени;
# [1:4,5:5:20] - выполняется аппроксимация полиномами до 20-ой степени, и выводится результат после аппроксимации полиномами степени 1,2,3,4,5,10,15,20. Обозначения 1:4 и 5:5:20 имеют тот же смысл, что и для циклов в языке Julia.

# <windowsparams> - параметры окна обработки данных.
# В скрипте заложена возможность обработки данных тремя способами:
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
# **********************************************************
# Примеры указания параметров окна:
# 1:end - формируется одно окно размером в весь входной файл;
# 1:1000 - формируется одно окно размером в 1000 точек, начальная точка 
# совпадает с первой точкой данных, конечная совпадает с 1000-й точкой данных;
# 1:100:1000 - формируются окна размером 100 точек по следующим отрезкам: (1, 100), (100, 199), (199, 298), ..., (901, 1000);
# 1:100:10:1000 - формируются окна размером 100 точек c перекрытием 10 точек по следующим отрезкам: (1, 100), (91, 190), (181, 280), ..., (901, 1000).

############################################################
# РЕЗУЛЬТАТЫ РАБОТЫ СКРИПТА
############################################################

# В результате работы скрипта в директории запуска формируется поддиректория APPROX, 
# в которой создается еще одна поддиректория с именем <out_dir_name>,
# в <out_dir_name> по отдельности для каждого аппроксимируемого столбца из <read_columns>
# и каждой степени из <output_degrees> и каждого окна обработки 
# создается файл из 4 столбцов: 
# <дискретный аргумент>  
# <исходные аппроксимируемые значения>  
# <значения аппроксимирующего полинома>  
# <невязка аппроксимации = разность второго и третьего стобца>
# В имени файла указываются аппроксимируемый столбец, начало окна аппроксимации (номер строки файла,
# не считая комментариев и пропущенных строк) и степень аппроксимирующего полинома.
# Все числа внутри файла выводятся @sprintf("%1.15e",...) - можно увеличить точность при необходимости.
# Для удобства дальнейшего построения графиков с помощью gnuplot расширение всех файлов: *.gplot_4
# **********************************************************
# Дополнительно в директории APPROX формируется файл result_max_resid.txt  (путем дописывания новых строк в конец, если этот файл уже есть),
# в котором содержится краткая информация о результате аппроксимации:
# максимум невязки при максимальной степени аппроксимации окна для каждого из входных столбцов
# В строку записываются следующие данные:
# - <out_dir_name> имя выходной директории аппроксимации
# - <current_window> обрабатываемое окно
# - максимум невязки для каждого аппроксимируемого столбца на текущем окне при максимальной степени аппроксимации

# **********************************************************
# Для логирования ошибок, если они возникли в ходе работы скрипта, в директории запуска скрипта создается файл errorlog_approx.log,
# в который выводятся следующие виды ошибок с указанием названия вычислительного эксперимента и времени запуска:
# "Некорректные данные, считано 0 строк" - происходит в случае, если скрипт не смог 
# прочитать исходный файл или были указаны неверные номера столбцов;
# "Некорректные данные, длина окна больше, чем размер данных" - 
# информирует о том, что при задании параметров окна был неверно выбран размер окна 
# и данных недостаточно, чтобы покрыть хотя бы 1 окно.

############################################################
# Пример запуска: 
############################################################
# julia.exe <script.jl> F-testdata.txt [1:4] [0,1,3,5:5:20] 1:end
# где:
# 5:5:20 в <output_degrees> означает список степеней с началом 5, шагом 5 и концом 20
# 1:end означает обработку всего интервала "меток времени", с начала файла до конца как одного окна
############################################################

############################################################
############################################################

if length(ARGS) != 4
    println(length(ARGS))
    println("warning: Неверное количество входных аргументов")
    exit()
end

using DelimitedFiles
using LinearAlgebra
using Printf
using Dates

# Подключение модуля вспомогательных функций

include("service_functions-ru.jl")

# Начало основного скрипта аппроксимации

out_dir_name = "equ_$(ARGS[1])_col$(ARGS[2])_degs$(ARGS[3])_w[$(ARGS[4])]_Float64"
out_dir_name = replace(out_dir_name, ":"=>"-")

errorlog = "errorlog_approx.log"
max_resid_file = "APPROX/result_max_resid.txt"

filename = ARGS[1]
read_columns = make_numberlist(ARGS[2])

DATA = readbigfile(filename, read_columns, Float64)
length_DATA = length(DATA[:,1])

if length_DATA == 0
    add_log_messange(errorlog, "$(out_dir_name) Некорректные данные, считано 0 строк")
    exit()
end

# Проверка на равномерность шага "меток времени":
# определим значение первого аргумента
pre_value_argument = DATA[1,1]
# определим шаг в аргументах
step_value_argument = DATA[2,1] - DATA[1,1]
# точность операции сравнения
prec_toler = step_value_argument * 1e-10

for ind_data = 2:length_DATA
    global pre_value_argument
    pva=pre_value_argument
    if abs((DATA[ind_data,1]-pva) - step_value_argument)>prec_toler
        println("warning: Данные не эквидистантны, используйте скрипт для не эквидистантных данных")
        exit()
    end
    pre_value_argument = DATA[ind_data,1]
end

output_degrees = make_numberlist(ARGS[3])
windows_list = make_windows(ARGS[4], length_DATA)

if length(windows_list) == 0
    add_log_messange(errorlog, "$(out_dir_name) Некорректные данные, длина окна больше, чем размер данных")
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

# Подготовка таблицы значений дискретных ортогональных полиномов Хана необходимых степеней на основе первого окна
# в массиве pnts - значения дискретного аргумента  (задаваемого первой прочитанной колонкой) одномерной дискретной (вектор)-функции (задаваемой остальными прочитанными колонками)
pnts = DATA[windows_list[1][1]:windows_list[1][2],1]
NN = length(pnts)
p1 = pnts[1]
pNN = pnts[end]
# нормируем дискретный аргумент на интервал [-1,1]:
pnts11 = zeros(Float64, NN)
for k=1:NN
    pnts11[k] = 2*(pnts[k] - p1) / (pNN - p1) - 1
end
# pnts11 = 2 * (pnts) / (pNN - p1) - 1
# Создание полиномов
pols = zeros(Float64, NN, maxdeg + 1)
# в массиве pols[j,k] будет формироваться значение ортог. полинома степени (k-1) в j-й точке дискретного аргумента (задаваемого первой прочитанной колонкой)
# Вначале в этот массив будут внесены значения ортогональных полиномов Лежандра (они НЕ ортогональны на нашей дискретной решетке аргумента!), а потом специальным способом ортогонализованы и нормированы (вариант алгоритма Грама-Шмидта).
# Первый полином - константа
pols[:,1] = ones(Float64, NN)
# второй - линейная функция дискретного аргумента (задаваемого первой прочитанной колонкой)
pols[:,2] = pnts11
# далее формируем (не дискретно-ортогональные!) ряды значений полиномов Лежандра (не нормированных) в точках pnts11
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
# полученные значения полиномов Лежандра теперь надо ортогонализовать!
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
            # - именно такой вариант оказался устойчив при ортогонализации, даже для Float64 !!!!
        end
        nev = 0.0
        for i = 1:(n - 1)
            nev = max(nev, abs(sum(u .* pols[:,i])))
        end
    end
    norm0 = sqrt(sum(abs2, u))
    pols[:,n] = u / norm0
    println("degree = ", n - 1, "   max nev after orthogonalization = ", @sprintf("%1.15e",nev),
    "   repN = ", repN, "   norm was: ", @sprintf("%1.15e",norm0))
    flush(stdout)
end
# Конец ортогонализации - в массиве pols[:,:] имеем набор значений дискретных ортогональных полиномов Хана,
#  ортогональных и нормированных на множестве точек pnts
println("Orthogonalization done")
flush(stdout)


# Обработка окон в файле

for current_window in windows_list
    println("Processing window: $current_window")
    flush(stdout)
    max_resid_in_col = zeros(Float64, length(read_columns)-1)
    for ind_cn = 2:length(read_columns)
        data_col = DATA[current_window[1]:current_window[2],ind_cn]
        data_arg = DATA[current_window[1]:current_window[2], 1]
        data_cn = copy(data_col)
        pol_approx = zeros(Float64, NN)
        # аппроксимация данных (по отдельности каждого столбца) полиномами - поскольку они ортогональны и нормированы, достаточно простого скалярного произведения!
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
                close(out_file)
                println(" written!")
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
