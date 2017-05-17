rm result.gal
#time ./galsim 10 input_data/ellipse_N_00010.gal 200 1e-5 0
#./test 10 result.gal ref_output_data/ellipse_N_00010_after200steps.gal
# time ./galsim 100 input_data/ellipse_N_00100.gal 200 1e-5 0
# ./test 100 result.gal ref_output_data/ellipse_N_00100_after200steps.gal
 time ./galsim 500 input_data/ellipse_N_00500.gal 200 1e-5 0
 ./test 500 result.gal ref_output_data/ellipse_N_00500_after200steps.gal
# time ./galsim 1000 input_data/ellipse_N_01000.gal 200 1e-5 0
# ./test 1000 result.gal ref_output_data/ellipse_N_01000_after200steps.gal
# time ./galsim 2000 input_data/ellipse_N_02000.gal 200 1e-5 0
# ./test 2000 result.gal ref_output_data/ellipse_N_02000_after200steps.gal



# ./test 10 result.gal ref_output_data/ellipse_N_00010_after200steps.gal
# ./test 100 result.gal ref_output_data/ellipse_N_00100_after200steps.gal
# ./test 500 result.gal ref_output_data/ellipse_N_00500_after200steps.gal
# ./test 1000 result.gal ref_output_data/ellipse_N_01000_after200steps.gal
# ./test 2000 result.gal ref_output_data/ellipse_N_02000_after200steps.gal
