#!/bin/bash

num_of_iterations=400000

cd k_0_1282
./k_0_1282.out 0.00 20.0 34.0 $num_of_iterations
cd ..
cd k_0_2564
./k_0_2564.out 0.00 20.0 34.0 $num_of_iterations
cd ..
cd k_0_641
./k_0_641.out 0.00 20.0 34.0 $num_of_iterations
cd ..
cd k_1_4102
./k_1_4102.out 0.00 20.0 34.0 $num_of_iterations
cd ..
cd k_1_923
./k_1_923.out 0.00 20.0 34.0 $num_of_iterations
cd ..
cd k_2_4358
./k_2_4358.out 0.00 20.0 34.0 $num_of_iterations
cd ..
cd k_first_two
./k_first_two.out 0.00 20.0 34.0 $num_of_iterations
cd ..
cd k_first_five
./k_first_five.out 0.00 20.0 34.0 $num_of_iterations
cd ..
cd k_first_eleven
./k_first_eleven.out 0.00 20.0 34.0 $num_of_iterations
cd ..
cd k_first_fifteen
./k_first_fifteen.out 0.00 20.0 34.0 $num_of_iterations
cd ..
cd k_first_nineteen
./k_first_nineteen.out 0.00 20.0 34.0 $num_of_iterations
cd ..
