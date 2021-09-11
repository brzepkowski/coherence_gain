# For measuring time gains
python3 -m timeit -s "from coherence_gain_t_finite_detailed import calc_W_one_part" "print('{:,}'.format(calc_W_one_part(t=0.2, T=72)))"
