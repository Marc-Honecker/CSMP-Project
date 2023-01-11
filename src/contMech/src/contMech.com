rm contMech.x
g++ -O2 *.cpp -std=c++11 -lfftw3 -L  ~/bin/fftw/lib -I  ~/bin/fftw/include -o contMech.x >& contMech.e
