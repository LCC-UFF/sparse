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

# create Coo_matrix
function coo_matrix(i::Vector{Int},j::Vector{Int},v::Vector)
    m = maximum(i)
    n = maximum(j)
    Coo_matrix(i,j,v,m,n)
end

# create Coo_matrix from existing dense matrix
function coo_matrix(mat::Matrix)
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

# add new value to existing Coo_matrix
function new_value!(mat::Coo_matrix,new::Tuple{Int,Int,Any})
    mat.m = max(mat.m,new[1])
    mat.n = max(mat.n,new[2])
    push!(mat.i,new[1])
    push!(mat.j,new[2])
    push!(mat.v,new[3])
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



################  CSR   ################

####### main struct definition
mutable struct Csr_matrix{T}
    row_pts::Vector{Int}
    collums::Vector{Int}
    values::Vector{T}
    m::Int
    n::Int    
end



####### Csr_matrix creation functions

# create Csr_matrix
function csr_matrix(i::Vector{Int},j::Vector{Int},v::Vector,m::Int=0,n::Int=0)
    m = max(maximum(i),m)
    n = max(maximum(j),n)
    row_pts = ones(Int,m+1)

    collums = [j[1]]
    values = [v[1]]
    row_pts[i[1]+1:m+1] = row_pts[i[1]+1:m+1] .+ 1

    for index in 2:size(v,1)
        insert!(collums,row_pts[i[index]+1],j[index])
        insert!(values,row_pts[i[index]+1],v[index])
        row_pts[i[index]+1:m+1] = row_pts[i[index]+1:m+1] .+ 1
    end
    Csr_matrix(row_pts,collums,values,m,n)
end

# create Csr_matrix from existing dense matrix
function csr_matrix(mat::Matrix)
    m = size(mat,1)
    n = size(mat,2)

    row_pts = ones(Int,m+1)
    collums = Vector{Int}(undef,0)
    values = Vector{eltype(mat)}(undef,0)
    
    for i in 1:m
        count = 0
        for j in 1:n
            if mat[i,j] != 0
                push!(collums,j)
                push!(values,mat[i,j])
                count += 1
            end
        end
        row_pts[i+1:m+1] = row_pts[i+1:m+1] .+ count
    end
    Csr_matrix(row_pts,collums,values,m,n)
end

# add new value to existing Csr_matrix
function new_value!(mat::Csr_matrix,new::Tuple{Int,Int,Any})
    mat.m = max(new[1],mat.m)
    mat.n = max(new[2],mat.n)
    insert!(mat.collums,mat.row_pts[new[1]+1],new[2])
    insert!(mat.values,mat.row_pts[new[1]+1],new[3])
    mat.row_pts[new[1]+1:mat.m+1] = mat.row_pts[new[1]+1:mat.m+1] .+ 1
end


####### Related functions

# read element value
function getvalue(mat::Csr_matrix,x::Int,y::Int)
    value = 0
    for i in mat.row_pts[x]:mat.row_pts[x+1]-1
        if y == mat.collums[i]
            value += mat.values[i]
        end
    end
    return value
end
# read element value as dense matrix (a[i,j])
Base.getindex(mat::Csr_matrix, x::Int, y::Int) = getvalue(mat,x,y)

# calculates a matrix-vector product
function matvecprod(mat::Csr_matrix,vec::Vector)
    result_vec = zeros(mat.m)

    for i in 1:mat.m
        for j in mat.row_pts[i]:mat.row_pts[i+1]-1
            result_vec[i] += mat.values[j]*vec[mat.collums[j]]
        end
    end

    return result_vec
end
# calculates a matrix-vector product using only x*y
Base.:*(x::Csr_matrix, y::Vector) = matvecprod(x,y)