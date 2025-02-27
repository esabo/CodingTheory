using CodingTheory
using Oscar
using LinearAlgebra
using BenchmarkTools

function visit_all_subsets(mutate, g_len, v_weight, subset)
  if subset
    gi=CodingTheory.SubsetGrayCode(g_len, v_weight)
  else
    gi=CodingTheory.GrayCode(g_len, v_weight; mutate = mutate)
  end
  for subset in gi
    # BenchmarkTools can elide simple computations so we need to do some nontrivial calculation here.
    # In any realistic situation we need to look at least look at all entries of the vector. 
    # Ill model the task of looking at the entries using the 'all' function
    Base.all(i->(i==0), subset) 
  end
  return 
end

#=
SubsetGrayCode and GrayCode appear to be in the same ballpark. 
GrayCode is faster, SubsetGrayCode uses less memory.

######
# n=25, k=12, mutate=true
######
@benchmark visit_all_subsets(25, 12, false)
BenchmarkTools.Trial: 54 samples with 1 evaluation.
 Range (min … max):  66.365 ms … 219.557 ms  ┊ GC (min … max):  0.00% … 34.15%
 Time  (median):     92.912 ms               ┊ GC (median):    26.24%
 Time  (mean ± σ):   93.526 ms ±  27.981 ms  ┊ GC (mean ± σ):  19.78% ± 13.74%

  ▂▂            █▃                                              
  ███▅▅▁▁▁▁▄▁▁█████▁▁▄▁▁▁▁▄▁▁▁▁▁▅▁▁▁▁▁▁▁▁▁▄▁▁▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄ ▁
  66.4 ms         Histogram: frequency by time          181 ms <

 Memory estimate: 158.70 MiB, allocs estimate: 5200322.

BenchmarkTools.Trial: 65 samples with 1 evaluation.
 Range (min … max):  51.159 ms … 103.931 ms  ┊ GC (min … max):  0.00% … 50.88%
 Time  (median):     77.176 ms               ┊ GC (median):    31.98%
 Time  (mean ± σ):   77.271 ms ±   6.824 ms  ┊ GC (mean ± σ):  31.52% ±  7.35%

                                          █ ▁                   
  ▃▃▁▁▁▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▅▆████▅▆▅▁▃▃▁▁▄▁▁▃▁▁▁▁▃ ▁
  51.2 ms         Histogram: frequency by time         89.3 ms <

 Memory estimate: 238.05 MiB, allocs estimate: 5200308.
=#