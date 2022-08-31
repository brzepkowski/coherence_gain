#!/bin/bash

num_of_iterations=1000

cd k=0.1282
./k=0.1282.out 0.00 20.0 34.0 $num_of_iterations
cd ..
cd k=0.2564
./k=0.2564.out 0.00 20.0 34.0 $num_of_iterations
cd ..
cd k=0.641
./k=0.641.out 0.00 20.0 34.0 $num_of_iterations
cd ..
cd k=1.4102
./k=1.4102.out 0.00 20.0 34.0 $num_of_iterations
cd ..
cd k=1.923
./k=1.923.out 0.00 20.0 34.0 $num_of_iterations
cd ..
cd k=2.4358
./k=2.4358.out 0.00 20.0 34.0 $num_of_iterations
cd ..
cd k=first_two
./k=first_two.out 0.00 20.0 34.0 $num_of_iterations
cd ..
cd k=first_five
./k=first_five.out 0.00 20.0 34.0 $num_of_iterations
cd ..
cd k=first_eleven
./k=first_eleven.out 0.00 20.0 34.0 $num_of_iterations
cd ..
cd k=first_fifteen
./k=first_fifteen.out 0.00 20.0 34.0 $num_of_iterations
cd ..
cd k=first_nineteen
./k=first_nineteen.out 0.00 20.0 34.0 $num_of_iterations
cd ..
