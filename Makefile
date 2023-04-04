all: legendre_otc cont_frac

legendre_otc: test_legendre_otc.cpp legendre_otc.h
	g++ test_legendre_otc.cpp -o test_legendre_otc -lstdc++

gsl_test: test_gsl.cpp
	g++ test_gsl.cpp -o test_gsl -lgsl -lgslcblas

cont_frac: test_cont_frac.cpp legendre_otc.h
	g++ test_cont_frac.cpp -o test_cont_frac -lstdc++
