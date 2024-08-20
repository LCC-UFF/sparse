################  COO   ################

####### main struct definition
mutable struct Coo_matrix{T}
    i::Vector{Int}
    j::Vector{Int}
    v::Vector{T}
    m::Int
    n::Int
end



####### Coo_matrix creation functions

# create Coo_matrix without defined matrix size
function Coo_matrix(i::Vector,j::Vector{Int},v::Vector)
    m = maximum(i)
    n = maximum(j)
    Coo_matrix(i,j,v,m,n)
end

# create Coo_matrix from existing dense matrix
function Coo_matrix(mat::Matrix)
    i_vec::Vector{Int} = []
    j_vec::Vector{Int} = []
    v_vec::Vector{eltype(mat)} = []
    for i in 1:size(mat,1)
        for j in 1:size(mat,2)
            if mat[i,j] != 0
                push!(i_vec,i)
                push!(j_vec,j)
                push!(v_vec,mat[i,j])
            end
        end
    end
    return Coo_matrix(i_vec,j_vec,v_vec,size(mat,1),size(mat,2))
end



####### Related functions

# read element value
function getvalue(mat::Coo_matrix,x::Int,y::Int)
    value = 0
    for i in 1:size(mat.i,1)
        if x==mat.i[i]
            if y==mat.j[i]
                value += mat.v[i]
            end
        end
    end
    return value
end
# read element value as dense matrix (a[i,j])
Base.getindex(mat::Coo_matrix, x::Int, y::Int) = getvalue(mat,x,y)

# calculates a matrix-vector product
function matvecprod(mat::Coo_matrix,vec::Vector)
    result_vec = zeros(mat.m)
    for count in 1:size(mat.v,1)
        result_vec[mat.i[count]] += mat.v[count]*vec[mat.j[count]]
    end
    return result_vec
end
# calculates a matrix-vector product using only x*y
Base.:*(x::Coo_matrix, y::Vector) = matvecprod(x,y)
