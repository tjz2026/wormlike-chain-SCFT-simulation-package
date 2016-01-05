rm *.o test
mpicxx -c  -g test.cpp global.cpp fft.cpp initial.cpp G_matrix.cpp MDE_linear.cpp cell.cpp field.cpp density_energy.cpp scft_io.cpp AB_diblock_driver.cpp matrix.cpp m_f.cpp cinv.cpp -I /home/tangjz/Soft/fftw3.3.3-double/include/
mpicxx *.o -o test -L /home/tangjz/Soft/fftw3.3.3-double/lib/ -lfftw3_mpi -lfftw3 -lm
