About
-----

A Go wrapper for the Levmar library which implements the Levenberg-Marquardt non-linear least squares algorithm

www.ics.forth.gr/~lourakis/levmar/

depends on 
- C libraries:   lapack, blas, f2c
- Go libraries:  github.com/verdverm/{go-symexpr,go-pge/problems}


Building
--------

- sudo apt-get install liblapack3 liblapack-dev libblas3 libblas-dev f2c 
- cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DLINSOLVERS_RETAIN_MEMORY=0 .
- make

