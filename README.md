# High degree least squares polynomial fitting using discrete orthogonal polynomials

© Siberian Federal University, 2017

## Legal disclaimer

Dear user!

Siberian Federal University recognizes the value of science and is committed to integration into the global world and dissemination of knowledge.

Guided by the principles of the open scientific community, SibFU decided to post the most complete information about some of its achievements for the development of technologies around the world.

At the same time, following the standards of international law and Russian law, the Siberian Federal University recognizes, secures and protects its own exclusive rights and the right of authorship of our researchers to the created intellectual results.

In connection with the foregoing, before using our programs, we ask you to read the terms of a free (non-exclusive) open license in the file `SFU OPEN LICENSE 2020.docx` for our products and join it.

Being a state institution, the Siberian Federal University has the obligation to properly account for programs and the facts of their use. Please kindly fill out a simple consent (acceptance) confirming your use of the program. If the legislation of the territory of your research binds for this, certify the signature in your scientific organization and send a scanned copy of the consent to the email address ykrot@sfu-kras.ru.

## Getting Started

This repository contains the code used in:
S.P.Tsarev, A.A.Kytmanov, _Discrete orthogonal polynomials as a tool for detection of small anomalies of time series: a case study of GPS final orbits_ (Posted at: <https://arxiv.org/> in 2020)

All algorithms for best least squares polynomial approximation in this repository were implemented in [Julia programming language](https://en.wikipedia.org/wiki/Julia_(programming_language)). 

The Julia scripts robustly approximate a (probably vector-valued) time series with a very large number of points (thousands) with polynomials of large degrees (up to a few hundred) 
in the sense of the best RMS (best least squares) approximation
and find the approximation residue.
You can approximate your dataset by polynomials on the entire data (or a part thereof) available in the file,
or alternatively approximate the data using a sliding window (possibly with overlapping).

Details on script parameters and the obtained results are given as comments inside the scripts themselves.

The scripts come in 4 different versions:
- for calculation with the standard (8-byte) floating point type `Float64` for time series with constant time step (equidistant argument lattice);
- for calculation with the standard (8-byte) floating point type `Float64` for time series with (slightly) non-constant time step (nonequidistant argument lattice);
- for calculation with the arbitrary precision floating point type `BigFloat` for time series with constant time step (equidistant argument lattice);
- for calculation with the arbitrary precision floating point type `BigFloat` for time series with (slightly) non-constant time step (nonequidistant argument lattice);

When you are not sure if your problem may involve fast error accumulation use `BigFloat` versions, but they are _slow_ and more memory-consuming!
Usually (for data up to 10,000 points) and approximation up to degree 200 you will be safe with `Float64` versions.
We found that using a modest size sliding window on the big time series may be more informative and much faster than approximation of the complete set of data: 
the log file `result_max_resid.txt` will give a concise report on approximation results for every window.

The Julia scripts for processing time series data (with English comments) can be found in subdirectory [**en**](en).

Complete examples of running the scripts with the log files and results can be found in subdirectory [**test**](test).

The Julia scripts with comments in Russian (otherwise absolutely identical to the English versions) as well as the `Legal disclaimer (ru).txt` and `SFU OPEN LICENSE 2020 (ru).docx` from Siberian Federal University can be found in subdirectory [**ru**](ru).

### Installing

You need only a Julia installation to run the code. Please note that all version require the file `service_functions-en.jl` (or `service_functions-ru.jl` for Russian language versions) where some utility functions are collected.

## Authors

[Sergey P. Tsarev](https://scholar.google.ru/citations?hl=ru&user=dey5ZnMAAAAJ&view_op=list_works), Alexander S. Pustoshilov

Please write to <sptsarev@mail.ru> if you have questions, comments or suggestions.

Cite the paper 
S.P.Tsarev, A.A.Kytmanov, _Discrete orthogonal polynomials as a tool for detection of small anomalies of time series: a case study of GPS final orbits_ 
if you need a publication reference.
