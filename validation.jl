include("spvec.jl")


#### teste coo_matrix

# a = coo_matrix([1,1],[1,2],[1,4])

# new_value!(a,(2,1,3))

# println(a)
# println(a[1,1])

# b = [1 1 0;
#      2 0 2]

# c = coo_matrix(b)
# println(c)
# println(c[2,1])

# d = [2,2,2]

# println(c*d)



#### teste csr_matrix

a = csr_matrix([1,2,1,3],[1,2,3,1],[1,2,3,4])

new_value!(a,(2,1,3))

println(a)
println(a[3,1])

b = [1 1 0;
     2 0 2]

c = csr_matrix(b)
println(c)
println(c[2,1])

d = [2,2,2]

println(c*d)