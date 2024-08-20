include("spvec.jl")


a = Coo_matrix([1,1],[1,2],[1,4])

println(a)
println(a[1,1])

b = [1 1 0;
     2 0 2]

c = Coo_matrix(b)
println(c)
println(c[2,1])

d = [2,2,2]

println(c*d)
