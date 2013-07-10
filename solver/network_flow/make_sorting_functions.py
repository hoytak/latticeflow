from math import ceil, log, floor
import random

output_file   = "fixed_sorting_functions.hpp"
top_network_n = 128

# file_template needs to have %(functions)s
# full_template needs to have %(swaps)s and %(size)d
# swap_template needs to have %(idx1)d. and %(idx2)d.


file_template = open('sort_function_templates/sorting_functions_file.cpt').read()
full_template = open('sort_function_templates/sorting_functions_func.cpt').read()
swap_template = open('sort_function_templates/sorting_functions_swap.cpt').read()

def generate_sorting_list(n):

    swap_instructions = []


    def log2_lt(n):
    
        k=1
    
        while k<n:
            k *= 2
            
        return k / 2
    
    def bitonicSort(start_idx, n, up):
        
        if n > 1:
        
            m = int(n/2)
            
            bitonicSort(start_idx, m, not up)
            
            bitonicSort(start_idx+m, n-m, up)
            
            bitonicMerge(start_idx, n, up)

    def bitonicMerge(start_idx, n, up):
        if n > 1:
            m = log2_lt(n)
            
            for i in range(start_idx, start_idx + n - m):
                swap_instructions.append( (i, i + m) if up else (i + m, i) )        
                
            bitonicMerge(start_idx, m, up)
            bitonicMerge(start_idx+m, n-m, up)

    bitonicSort(0, n, True)

    return swap_instructions


def test_sn(n):
    
    sl = generate_sorting_list(n)

    for iteration in range(1000):
        x = [random.random() for i in range(n)]
        
        for i1, i2 in sl:
            if x[i1] > x[i2]:
                x[i1], x[i2] = x[i2], x[i1]
            
        if not x == sorted(x):
            print "ERROR!!!" 
            print "nsx = ", x
            print "xl = ", sorted(x)
            print "mask = ", ''.join('-' if a == b else 'X' for a,b in zip(x, sorted(x)))
            print "diff = ", [a - b for a,b in zip(x, sorted(x))]
            print sl

            return False
    print "Testing for n = %d passes." % n
    return True
    
    
#for i in range(2,500):
#    if not test_sn(i):
#        break

def getSortingCode(n):
    """
    Gets deterministic network based sorting code for small arrays.
    
    """
    
    sl = generate_sorting_list(n)
        
        
    return 
if __name__ == '__main__':
    
    full_file = file_template % {"functions" 
                                 : "\n\n".join(
                                     (full_template 
                                      % {'size' : n,
                                         'swaps' : '\n'.join((swap_template 
                                                              % { 'idx1' : idx1,
                                                                  'idx2' : idx2 })
                                                             for idx1, idx2 in generate_sorting_list(n))})
                                     for n in range(2, top_network_n))} 
        
    
    f = open(output_file, 'w')
    f.write(full_file)
    f.close()
    
    print "sorting functions written to file %s." % output_file
    
    

