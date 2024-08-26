# Compute Effective ELASTIC Properties
using JSON
using SparseArrays
using LinearAlgebra

# Model data struct:
struct Model
    nx::UInt64;
    ny::UInt64;
    nz::UInt64;
    voxelSize::Float64;
    refinement::UInt64;
    nMat::UInt16;
    rhsType::UInt8;
    solverType::UInt8;
    pcgTol::Float64;
    pcgIter::UInt64;
    matKeys::Vector{UInt16};
    matProp::Matrix{Float64}; 
    nNodes::UInt64;
    nElems::UInt64;
    nDOFs::UInt64;
    DOFMap::Vector{UInt64}; 
    elemMatMap::Vector{UInt16}; 
    function Model(_nx, _ny, _nz, _voxelSize, _refinement, _nMat, _rhsType, _solverType, _pcgTol, _pcgIter, _matKeys, _matProp, _nNodes, _nElems, _nDOFs, _DOFMap, _elemMatMap)
        new(_nx, _ny, _nz, _voxelSize, _refinement, _nMat, _rhsType, _solverType, _pcgTol, _pcgIter, _matKeys, _matProp, _nNodes, _nElems, _nDOFs, _DOFMap, _elemMatMap);
    end
end

# Build Model:
function buildModel(_JsonFile::String, _RawFile::String)
    println(".Building Model!")
    # Read Json file:
    nx, ny, nz, voxelSize, refinement, nMat, rhsType, solverType, pcgTol, pcgIter, matKeys, matProp = readJSON(_JsonFile)
    # Read Raw file:
    elemMatMap = zeros(UInt16, nx * ny * nz * refinement * refinement * refinement)
    readRAW!(nx, ny, nz, refinement, matKeys, elemMatMap, _RawFile)
    # Update the parameters based on the given refinement level:
    nx *= refinement
    ny *= refinement
    nz *= refinement
    nNodes::UInt64 = (nx + 1) * (ny + 1) * (nz + 1)
    nElems::UInt64 = (nx) * (ny) * (nz)
    DOFperNode::UInt64 = 3
    nDOFs::UInt64 = nElems * DOFperNode
    # Generate a map of Degree of Freedom:
    DOFMap = zeros(UInt64, nNodes)
    generateDOFMap!(nx, ny, nz, DOFMap)   
    # Build the Model:
    model = Model(nx, ny, nz, voxelSize, refinement, nMat, rhsType, solverType, pcgTol, pcgIter, matKeys, matProp, nNodes, nElems, nDOFs, DOFMap, elemMatMap)
    println("---------------------------")
    return model
end

# Read JSON file:
function readJSON(_filename::String)
    println("   .Reading JSON!")
    # Open and read file:
    open(_filename, "r") do f
        data = JSON.parse(f)
        nx::UInt64 = data["image_dimensions"][1]
        ny::UInt64 = data["image_dimensions"][2]
        nz::UInt64 = data["image_dimensions"][3]
        refinement::UInt64 = 1
        if haskey(data, "refinement");       refinement = data["refinement"];     end
        voxelSize::Float64 = 1.0
        if haskey(data, "voxel_size"); voxelSize = data["voxel_size"]; end
        rhsType::UInt8 = 0
        if haskey(data, "type_of_rhs");      rhsType = data["type_of_rhs"];       end
        solverType::UInt8 = 0
        if haskey(data, "type_of_solver");   solverType = data["type_of_solver"]; end
        pcgTol::Float64 = 0.000001
        if haskey(data, "solver_tolerance"); pcgTol = data["solver_tolerance"];   end       
        pcgIter::UInt64 = nx * ny * nz * refinement * refinement * refinement
        if haskey(data, "number_of_iterations"); pcgIter = data["number_of_iterations"]; end
        nMat::UInt16 = data["number_of_materials"]
        materials = data["properties_of_materials"]
        matKeys = zeros(UInt16, 256)
        matProp = zeros(Float64, 256, 2)
        for i = 1:nMat
            matKeys[convert(UInt8, materials[i][1]) + 1] = i
            matProp[convert(UInt8, materials[i][1]) + 1,1] = convert(Float64, materials[i][2])
            matProp[convert(UInt8, materials[i][1]) + 1,2] = convert(Float64, materials[i][3])
        end
        materials = nothing
        data = nothing
        return nx, ny, nz, voxelSize, refinement, nMat, rhsType, solverType, pcgTol, pcgIter, matKeys, matProp
    end    
end

# Read RAW file:
function readRAW!(_nx::UInt64, _ny::UInt64, _nz::UInt64, _refinement::UInt64, _matKeys::Vector{UInt16}, _elemMatMap::Vector{UInt16}, _filename::String)
    println("   .Reading RAW!")
    # Initializations
    nelem::UInt64 = _nx * _ny * _nz; slice::UInt64 = _nx * _ny; buffer::UInt64 = 0;    
    elref::UInt64 = 0; ix::UInt64 = 0; iy::UInt64 = 0; iz::UInt64 = 0;
    row::UInt64 = _ny * _refinement; rowref::UInt64 = _ny * _refinement * _refinement; 
    sliceref2::UInt64 = _refinement*_refinement*slice; 
    sliceref3::UInt64 = _refinement*_refinement*_refinement*slice; 
    # Open and read file
    open(_filename, "r") do io
        bin_array = read(io)      
        # Build the element material map based on the refinement level:
        for e = 1:nelem
            buffer = _matKeys[bin_array[e] + 1]
            ix = ((e-1) % _nx)
            iz = ((e-1) ÷ slice)
            iy = (((e-1)-iz*slice) ÷ _nx)   
            # el = 1 + (ix * _ny) + iy + (iz * slice)
            elref = 1 + (ix * rowref) + (iy*_refinement) + (iz*sliceref3) 
            for k = 1:_refinement
                for i = 1:_refinement
                    for j = 1:_refinement
                        _elemMatMap[elref + (j-1) + (i-1)*row + (k-1)*sliceref2] = buffer
                    end
                end        
            end 
        end  
        bin_array = nothing
    end
end

# Generate the Degree of Freedom Map:
function generateDOFMap!(_nx::UInt64, _ny::UInt64, _nz::UInt64, _DOFMap::Vector{UInt64})
    println("   .Generating the Map of DOFs (Degrees of Freedom)!")
    # Number the DOFs following the nodes from top to bottom and left to right:
    nElemS::UInt64 = _nx * _ny
    nNodeS::UInt64 = (_nx + 1) * (_ny + 1)
    i::UInt64 = 0
    @fastmath @inbounds @simd for n = 1:(nNodeS * (_nz + 1))
        i = (n - 1) % nNodeS
        _DOFMap[n] = (i - (i ÷ (_ny + 1)) - _ny * ((i % (_ny + 1)) ÷ _ny)) % nElemS + (((n - 1) ÷ nNodeS) % _nz) * nElemS + 1
    end
end

# Estimate memory consuption:
function estimateMemory(_model::Model)
    println("   .Estimating memory!")
    # elemMatMap = 16 bits * nElems
    # DOFMap = 64 bits * nNodes
    # RHS = 64 bits * nDOFs
    # PCGM Solver   / _solver == 0 / M d x q = 4 * 64 bits * nDOFs
    # Direct Solver / _solver == 1 / K = 18 * 64 bits * nElems (rough sparse estimative)
    mem::Float64 = 0.0
    if (_model.solverType == 0)
        mem += (16 * _model.nElems + 64 * _model.nNodes + 5 * 64 * _model.nDOFs) / 8 / 1_000_000
    elseif (_model.solverType == 1)
        mem += (16 * _model.nElems + 64 * _model.nNodes + 2 * 64 * _model.nDOFs + 18 * 64 * _model.nElems) / 8 / 1_000_000
    end
    println("   $(_model.nDOFs) DOFs")
    println("   $mem MB")
    println("---------------------------")
end

# Compute the element stiffness matrix for each material:
function elementStiffnessMatrices!(_model::Model, _K::Array{Float64,3}, _B::Array{Float64,3})
    println("   .Computing each element stiffness matrix!")
    # Compute the matrices for each material:
    i::UInt64 = 0;
    for j = 1:256
        if (_model.matKeys[j] != 0)
            i += 1
            elemProps = _model.matProp[j,:]
            _K[:,:,i], _B[:,:,i] = C3D8ElementStiffness(elemProps)
        end
    end
end

# Element C3D8 Stiffness - FEM
function C3D8ElementStiffness(_elemProps::Vector{Float64})
    # Initializations
    k  = zeros(Float64, 24, 24)
    BC = zeros(Float64, 6, 24)
    C  = zeros(Float64, 6, 6)
    E = _elemProps[1]
    p = _elemProps[2]
    # Element coords
    x = [0.;1.;1.;0.;0.;1.;1.;0.]
    y = [0.;0.;1.;1.;0.;0.;1.;1.]
    z = [0.;0.;0.;0.;1.;1.;1.;1.]
    # Constitutive matrix
    E = E / ((1. + p) * (1. - 2 * p))
    C[1,1] = 1 - p;            C[1,2] = p;                C[1,3] = p;
    C[2,1] = p;                C[2,2] = 1 - p;            C[2,3] = p;
    C[3,1] = p;                C[3,2] = p;                C[3,3] = 1 - p;
    C[4,4] = (1 - 2 * p) / 2;  C[5,5] = (1 - 2 * p) / 2;  C[6,6] = (1 - 2 * p) / 2;
    C .*= E
    # Gauss Points and Weights
    gp = [-1.0 / sqrt(3) 1.0 / sqrt(3)]
    gw = [1.0 1.0]
    for i = 1:2
        r = gp[1,i]
        wx = gw[1,i]
        for j = 1:2
            s = gp[1,j]
            wy = gw[1,j]
            for l = 1:2
                t = gp[1,l]
                wz = gw[1,l]
                B, J = C3D8BMatrix(r, s, t, x, y, z)
                dJ = det(J)
                k  += B' * C * B * dJ * wx * wy * wz
                BC += C * B * dJ * wx * wy * wz
            end
        end
    end
    return k, BC
end

# C3D8BMatrix - FEM
function C3D8BMatrix(r::Float64, s::Float64, t::Float64, x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64})
    # Initializations
    B = zeros(Float64, 6, 24)
    X = [x'; y'; z']
    # Compute B matrix and Jacobian
    dN1dr = -(1 - s) * (1 - t) * .125; dN2dr =  (1 - s) * (1 - t) * .125; dN3dr =  (1 + s) * (1 - t) * .125; dN4dr = -(1 + s) * (1 - t) * .125;
    dN5dr = -(1 - s) * (1 + t) * .125; dN6dr =  (1 - s) * (1 + t) * .125; dN7dr =  (1 + s) * (1 + t) * .125; dN8dr = -(1 + s) * (1 + t) * .125;
    dN1ds = -(1 - r) * (1 - t) * .125; dN2ds = -(1 + r) * (1 - t) * .125; dN3ds =  (1 + r) * (1 - t) * .125; dN4ds =  (1 - r) * (1 - t) * .125;
    dN5ds = -(1 - r) * (1 + t) * .125; dN6ds = -(1 + r) * (1 + t) * .125; dN7ds =  (1 + r) * (1 + t) * .125; dN8ds =  (1 - r) * (1 + t) * .125;
    dN1dt = -(1 - r) * (1 - s) * .125; dN2dt = -(1 + r) * (1 - s) * .125; dN3dt = -(1 + r) * (1 + s) * .125; dN4dt = -(1 - r) * (1 + s) * .125;
    dN5dt =  (1 - r) * (1 - s) * .125; dN6dt =  (1 + r) * (1 - s) * .125; dN7dt =  (1 + r) * (1 + s) * .125; dN8dt =  (1 - r) * (1 + s) * .125;
    dN = [dN1dr dN2dr dN3dr dN4dr dN5dr dN6dr dN7dr dN8dr;
          dN1ds dN2ds dN3ds dN4ds dN5ds dN6ds dN7ds dN8ds;
          dN1dt dN2dt dN3dt dN4dt dN5dt dN6dt dN7dt dN8dt];
    J = dN * X'
    dNdx = J \ dN
    B[1,1]  = dNdx[1,1];  B[2,2]  = dNdx[2,1];  B[3,3]  = dNdx[3,1];
    B[1,4]  = dNdx[1,2];  B[2,5]  = dNdx[2,2];  B[3,6]  = dNdx[3,2];
    B[1,7]  = dNdx[1,3];  B[2,8]  = dNdx[2,3];  B[3,9]  = dNdx[3,3];
    B[1,10] = dNdx[1,4];  B[2,11] = dNdx[2,4];  B[3,12] = dNdx[3,4];
    B[1,13] = dNdx[1,5];  B[2,14] = dNdx[2,5];  B[3,15] = dNdx[3,5];
    B[1,16] = dNdx[1,6];  B[2,17] = dNdx[2,6];  B[3,18] = dNdx[3,6];
    B[1,19] = dNdx[1,7];  B[2,20] = dNdx[2,7];  B[3,21] = dNdx[3,7];
    B[1,22] = dNdx[1,8];  B[2,23] = dNdx[2,8];  B[3,24] = dNdx[3,8];
    B[4,1]  = dNdx[2,1];  B[5,1]  = dNdx[3,1];  B[6,2]  = dNdx[3,1];
    B[4,2]  = dNdx[1,1];  B[5,3]  = dNdx[1,1];  B[6,3]  = dNdx[2,1];
    B[4,4]  = dNdx[2,2];  B[5,4]  = dNdx[3,2];  B[6,5]  = dNdx[3,2];
    B[4,5]  = dNdx[1,2];  B[5,6]  = dNdx[1,2];  B[6,6]  = dNdx[2,2];
    B[4,7]  = dNdx[2,3];  B[5,7]  = dNdx[3,3];  B[6,8]  = dNdx[3,3];
    B[4,8]  = dNdx[1,3];  B[5,9]  = dNdx[1,3];  B[6,9]  = dNdx[2,3];
    B[4,10] = dNdx[2,4];  B[5,10] = dNdx[3,4];  B[6,11] = dNdx[3,4];
    B[4,11] = dNdx[1,4];  B[5,12] = dNdx[1,4];  B[6,12] = dNdx[2,4];
    B[4,13] = dNdx[2,5];  B[5,13] = dNdx[3,5];  B[6,14] = dNdx[3,5];
    B[4,14] = dNdx[1,5];  B[5,15] = dNdx[1,5];  B[6,15] = dNdx[2,5];
    B[4,16] = dNdx[2,6];  B[5,16] = dNdx[3,6];  B[6,17] = dNdx[3,6];
    B[4,17] = dNdx[1,6];  B[5,18] = dNdx[1,6];  B[6,18] = dNdx[2,6];
    B[4,19] = dNdx[2,7];  B[5,19] = dNdx[3,7];  B[6,20] = dNdx[3,7];
    B[4,20] = dNdx[1,7];  B[5,21] = dNdx[1,7];  B[6,21] = dNdx[2,7];
    B[4,22] = dNdx[2,8];  B[5,22] = dNdx[3,8];  B[6,23] = dNdx[3,8];
    B[4,23] = dNdx[1,8];  B[5,24] = dNdx[1,8];  B[6,24] = dNdx[2,8];
    return B, J
end

# Compute the RHS: Boundary or Domain, _rhsType: Boundary = 0 || Domain = 1
# _axis 0 = X || _axis 1 = Y || _axis 2 = Z || _axis 3 = XY || _axis 4 = XZ || _axis 5 = YZ ||
function computeRHS!(_model::Model, _RHS::Vector{Float64}, _axis::Int, _K::Array{Float64,3}, _B::Array{Float64,3})
    # Initializations
    pElemDOFNum = zeros(UInt64, 24)
    bound = zeros(UInt8, 4)
    N1::UInt64 = 0; N2::UInt64 = 0; N3::UInt64 = 0; N4::UInt64 = 0;
    N5::UInt64 = 0; N6::UInt64 = 0; N7::UInt64 = 0; N8::UInt64 = 0;
    nElemS::UInt64 = _model.nx * _model.ny
    nNodeS::UInt64 = (_model.nx + 1) * (_model.ny + 1)
    c::UInt64 = 0; r::UInt64 = 0;
    # Compute each RHS (_axis) based on boundary or domain data (_rhsType)
    if _model.rhsType == 1     # Boundary
        println("   .Computing RHS - Boundary!")
        e::UInt64 = 0
        delta::Float64 = 0.0
        if _axis == 0     # _axis 0 = X
            delta = _model.nx
            c = _model.nx
            bound[1] = 4; bound[2] = 7; bound[3] = 16; bound[4] = 19;
            for dz = 1:_model.nz
                for r = 1:_model.ny
                    e = r + (c - 1) * _model.ny + (dz - 1) * (nElemS);
                    N1 = 2 + (e - 1) % (nElemS) + div((e - 1) % (nElemS), _model.ny) + div(e - 1, nElemS) * nNodeS;
                    N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
                    N5 = N1 + nNodeS; N6 = N2 + nNodeS; N7 = N3 + nNodeS; N8 = N4 + nNodeS;
                    pElemDOFNum[1]  = _model.DOFMap[N1] * 3 - 2; pElemDOFNum[2]  = _model.DOFMap[N1] * 3 - 1; pElemDOFNum[3]  = _model.DOFMap[N1] * 3;
                    pElemDOFNum[4]  = _model.DOFMap[N2] * 3 - 2; pElemDOFNum[5]  = _model.DOFMap[N2] * 3 - 1; pElemDOFNum[6]  = _model.DOFMap[N2] * 3;
                    pElemDOFNum[7]  = _model.DOFMap[N3] * 3 - 2; pElemDOFNum[8]  = _model.DOFMap[N3] * 3 - 1; pElemDOFNum[9]  = _model.DOFMap[N3] * 3;
                    pElemDOFNum[10] = _model.DOFMap[N4] * 3 - 2; pElemDOFNum[11] = _model.DOFMap[N4] * 3 - 1; pElemDOFNum[12] = _model.DOFMap[N4] * 3;
                    pElemDOFNum[13] = _model.DOFMap[N5] * 3 - 2; pElemDOFNum[14] = _model.DOFMap[N5] * 3 - 1; pElemDOFNum[15] = _model.DOFMap[N5] * 3;
                    pElemDOFNum[16] = _model.DOFMap[N6] * 3 - 2; pElemDOFNum[17] = _model.DOFMap[N6] * 3 - 1; pElemDOFNum[18] = _model.DOFMap[N6] * 3;
                    pElemDOFNum[19] = _model.DOFMap[N7] * 3 - 2; pElemDOFNum[20] = _model.DOFMap[N7] * 3 - 1; pElemDOFNum[21] = _model.DOFMap[N7] * 3;
                    pElemDOFNum[22] = _model.DOFMap[N8] * 3 - 2; pElemDOFNum[23] = _model.DOFMap[N8] * 3 - 1; pElemDOFNum[24] = _model.DOFMap[N8] * 3;
                    for i = 1:24
                        for j in bound
                            _RHS[pElemDOFNum[i]] -= (_K[i,j,_model.elemMatMap[e]]) * delta
                        end
                    end
                end
            end
        elseif _axis == 1 # _axis 1 = Y
            delta = _model.ny
            r = 1
            bound[1] = 8; bound[2] = 11; bound[3] = 20; bound[4] = 23;
            for dz = 1:_model.nz
                for c = 1:_model.nx
                    e = r + (c - 1) * _model.ny + (dz - 1) * (nElemS);
                    N1 = 2 + (e - 1) % (nElemS) + div((e - 1) % (nElemS), _model.ny) + div(e - 1, nElemS) * nNodeS;
                    N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
                    N5 = N1 + nNodeS; N6 = N2 + nNodeS; N7 = N3 + nNodeS; N8 = N4 + nNodeS;
                    pElemDOFNum[1]  = _model.DOFMap[N1] * 3 - 2; pElemDOFNum[2]  = _model.DOFMap[N1] * 3 - 1; pElemDOFNum[3]  = _model.DOFMap[N1] * 3;
                    pElemDOFNum[4]  = _model.DOFMap[N2] * 3 - 2; pElemDOFNum[5]  = _model.DOFMap[N2] * 3 - 1; pElemDOFNum[6]  = _model.DOFMap[N2] * 3;
                    pElemDOFNum[7]  = _model.DOFMap[N3] * 3 - 2; pElemDOFNum[8]  = _model.DOFMap[N3] * 3 - 1; pElemDOFNum[9]  = _model.DOFMap[N3] * 3;
                    pElemDOFNum[10] = _model.DOFMap[N4] * 3 - 2; pElemDOFNum[11] = _model.DOFMap[N4] * 3 - 1; pElemDOFNum[12] = _model.DOFMap[N4] * 3;
                    pElemDOFNum[13] = _model.DOFMap[N5] * 3 - 2; pElemDOFNum[14] = _model.DOFMap[N5] * 3 - 1; pElemDOFNum[15] = _model.DOFMap[N5] * 3;
                    pElemDOFNum[16] = _model.DOFMap[N6] * 3 - 2; pElemDOFNum[17] = _model.DOFMap[N6] * 3 - 1; pElemDOFNum[18] = _model.DOFMap[N6] * 3;
                    pElemDOFNum[19] = _model.DOFMap[N7] * 3 - 2; pElemDOFNum[20] = _model.DOFMap[N7] * 3 - 1; pElemDOFNum[21] = _model.DOFMap[N7] * 3;
                    pElemDOFNum[22] = _model.DOFMap[N8] * 3 - 2; pElemDOFNum[23] = _model.DOFMap[N8] * 3 - 1; pElemDOFNum[24] = _model.DOFMap[N8] * 3;
                    for i = 1:24
                        for j in bound
                            _RHS[pElemDOFNum[i]] -= (_K[i,j,_model.elemMatMap[e]]) * delta
                        end
                    end
                end
            end
        elseif _axis == 2 # _axis 2 = Z
            delta = _model.nz
            dz = _model.nz
            bound[1] = 15; bound[2] = 18; bound[3] = 21; bound[4] = 24;
            for r = 1:_model.ny
                for c = 1:_model.nx
                    e = r + (c - 1) * _model.ny + (dz - 1) * (nElemS);
                    N1 = 2 + (e - 1) % (nElemS) + div((e - 1) % (nElemS), _model.ny) + div(e - 1, nElemS) * nNodeS;
                    N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
                    N5 = N1 + nNodeS; N6 = N2 + nNodeS; N7 = N3 + nNodeS; N8 = N4 + nNodeS;      
                    pElemDOFNum[1]  = _model.DOFMap[N1] * 3 - 2; pElemDOFNum[2]  = _model.DOFMap[N1] * 3 - 1; pElemDOFNum[3]  = _model.DOFMap[N1] * 3;
                    pElemDOFNum[4]  = _model.DOFMap[N2] * 3 - 2; pElemDOFNum[5]  = _model.DOFMap[N2] * 3 - 1; pElemDOFNum[6]  = _model.DOFMap[N2] * 3;
                    pElemDOFNum[7]  = _model.DOFMap[N3] * 3 - 2; pElemDOFNum[8]  = _model.DOFMap[N3] * 3 - 1; pElemDOFNum[9]  = _model.DOFMap[N3] * 3;
                    pElemDOFNum[10] = _model.DOFMap[N4] * 3 - 2; pElemDOFNum[11] = _model.DOFMap[N4] * 3 - 1; pElemDOFNum[12] = _model.DOFMap[N4] * 3;
                    pElemDOFNum[13] = _model.DOFMap[N5] * 3 - 2; pElemDOFNum[14] = _model.DOFMap[N5] * 3 - 1; pElemDOFNum[15] = _model.DOFMap[N5] * 3;
                    pElemDOFNum[16] = _model.DOFMap[N6] * 3 - 2; pElemDOFNum[17] = _model.DOFMap[N6] * 3 - 1; pElemDOFNum[18] = _model.DOFMap[N6] * 3;
                    pElemDOFNum[19] = _model.DOFMap[N7] * 3 - 2; pElemDOFNum[20] = _model.DOFMap[N7] * 3 - 1; pElemDOFNum[21] = _model.DOFMap[N7] * 3;
                    pElemDOFNum[22] = _model.DOFMap[N8] * 3 - 2; pElemDOFNum[23] = _model.DOFMap[N8] * 3 - 1; pElemDOFNum[24] = _model.DOFMap[N8] * 3;
                    for i = 1:24
                        for j in bound
                            _RHS[pElemDOFNum[i]] -= (_K[i,j,_model.elemMatMap[e]]) * delta
                        end
                    end
                end
            end
        elseif _axis == 3 # _axis 3 = XY
            delta = _model.ny
            r = 1
            bound[1] = 7; bound[2] = 10; bound[3] = 19; bound[4] = 22;
            for dz = 1:_model.nz
                for c = 1:_model.nx
                    e = r + (c - 1) * _model.ny + (dz - 1) * (nElemS);
                    N1 = 2 + (e - 1) % (nElemS) + div((e - 1) % (nElemS), _model.ny) + div(e - 1, nElemS) * nNodeS;
                    N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
                    N5 = N1 + nNodeS; N6 = N2 + nNodeS; N7 = N3 + nNodeS; N8 = N4 + nNodeS;
                    pElemDOFNum[1]  = _model.DOFMap[N1] * 3 - 2; pElemDOFNum[2]  = _model.DOFMap[N1] * 3 - 1; pElemDOFNum[3]  = _model.DOFMap[N1] * 3;
                    pElemDOFNum[4]  = _model.DOFMap[N2] * 3 - 2; pElemDOFNum[5]  = _model.DOFMap[N2] * 3 - 1; pElemDOFNum[6]  = _model.DOFMap[N2] * 3;
                    pElemDOFNum[7]  = _model.DOFMap[N3] * 3 - 2; pElemDOFNum[8]  = _model.DOFMap[N3] * 3 - 1; pElemDOFNum[9]  = _model.DOFMap[N3] * 3;
                    pElemDOFNum[10] = _model.DOFMap[N4] * 3 - 2; pElemDOFNum[11] = _model.DOFMap[N4] * 3 - 1; pElemDOFNum[12] = _model.DOFMap[N4] * 3;
                    pElemDOFNum[13] = _model.DOFMap[N5] * 3 - 2; pElemDOFNum[14] = _model.DOFMap[N5] * 3 - 1; pElemDOFNum[15] = _model.DOFMap[N5] * 3;
                    pElemDOFNum[16] = _model.DOFMap[N6] * 3 - 2; pElemDOFNum[17] = _model.DOFMap[N6] * 3 - 1; pElemDOFNum[18] = _model.DOFMap[N6] * 3;
                    pElemDOFNum[19] = _model.DOFMap[N7] * 3 - 2; pElemDOFNum[20] = _model.DOFMap[N7] * 3 - 1; pElemDOFNum[21] = _model.DOFMap[N7] * 3;
                    pElemDOFNum[22] = _model.DOFMap[N8] * 3 - 2; pElemDOFNum[23] = _model.DOFMap[N8] * 3 - 1; pElemDOFNum[24] = _model.DOFMap[N8] * 3;
                    for i = 1:24
                        for j in bound
                            _RHS[pElemDOFNum[i]] -= (_K[i,j,_model.elemMatMap[e]]) * delta
                        end
                    end
                end
            end
        elseif _axis == 4 # _axis 4 = XZ
            delta = _model.nz
            dz = _model.nz
            bound[1] = 13; bound[2] = 16; bound[3] = 19; bound[4] = 22;
            for r = 1:_model.ny
                for c = 1:_model.nx
                    e = r + (c - 1) * _model.ny + (dz - 1) * (nElemS);
                    N1 = 2 + (e - 1) % (nElemS) + div((e - 1) % (nElemS), _model.ny) + div(e - 1, nElemS) * nNodeS;
                    N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
                    N5 = N1 + nNodeS; N6 = N2 + nNodeS; N7 = N3 + nNodeS; N8 = N4 + nNodeS;
                    pElemDOFNum[1]  = _model.DOFMap[N1] * 3 - 2; pElemDOFNum[2]  = _model.DOFMap[N1] * 3 - 1; pElemDOFNum[3]  = _model.DOFMap[N1] * 3;
                    pElemDOFNum[4]  = _model.DOFMap[N2] * 3 - 2; pElemDOFNum[5]  = _model.DOFMap[N2] * 3 - 1; pElemDOFNum[6]  = _model.DOFMap[N2] * 3;
                    pElemDOFNum[7]  = _model.DOFMap[N3] * 3 - 2; pElemDOFNum[8]  = _model.DOFMap[N3] * 3 - 1; pElemDOFNum[9]  = _model.DOFMap[N3] * 3;
                    pElemDOFNum[10] = _model.DOFMap[N4] * 3 - 2; pElemDOFNum[11] = _model.DOFMap[N4] * 3 - 1; pElemDOFNum[12] = _model.DOFMap[N4] * 3;
                    pElemDOFNum[13] = _model.DOFMap[N5] * 3 - 2; pElemDOFNum[14] = _model.DOFMap[N5] * 3 - 1; pElemDOFNum[15] = _model.DOFMap[N5] * 3;
                    pElemDOFNum[16] = _model.DOFMap[N6] * 3 - 2; pElemDOFNum[17] = _model.DOFMap[N6] * 3 - 1; pElemDOFNum[18] = _model.DOFMap[N6] * 3;
                    pElemDOFNum[19] = _model.DOFMap[N7] * 3 - 2; pElemDOFNum[20] = _model.DOFMap[N7] * 3 - 1; pElemDOFNum[21] = _model.DOFMap[N7] * 3;
                    pElemDOFNum[22] = _model.DOFMap[N8] * 3 - 2; pElemDOFNum[23] = _model.DOFMap[N8] * 3 - 1; pElemDOFNum[24] = _model.DOFMap[N8] * 3;
                    for i = 1:24
                        for j in bound
                            _RHS[pElemDOFNum[i]] -= (_K[i,j,_model.elemMatMap[e]]) * delta
                        end
                    end
                end
            end
        elseif _axis == 5 # _axis 5 = YZ
            delta = _model.ny
            r = 1
            bound[1] = 9; bound[2] = 12; bound[3] = 21; bound[4] = 24;
            for dz = 1:_model.nz
                for c = 1:_model.nx
                    e = r + (c - 1) * _model.ny + (dz - 1) * (nElemS);
                    N1 = 2 + (e - 1) % (nElemS) + div((e - 1) % (nElemS), _model.ny) + div(e - 1, nElemS) * nNodeS;
                    N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
                    N5 = N1 + nNodeS; N6 = N2 + nNodeS; N7 = N3 + nNodeS; N8 = N4 + nNodeS;
                    pElemDOFNum[1]  = _model.DOFMap[N1] * 3 - 2; pElemDOFNum[2]  = _model.DOFMap[N1] * 3 - 1; pElemDOFNum[3]  = _model.DOFMap[N1] * 3;
                    pElemDOFNum[4]  = _model.DOFMap[N2] * 3 - 2; pElemDOFNum[5]  = _model.DOFMap[N2] * 3 - 1; pElemDOFNum[6]  = _model.DOFMap[N2] * 3;
                    pElemDOFNum[7]  = _model.DOFMap[N3] * 3 - 2; pElemDOFNum[8]  = _model.DOFMap[N3] * 3 - 1; pElemDOFNum[9]  = _model.DOFMap[N3] * 3;
                    pElemDOFNum[10] = _model.DOFMap[N4] * 3 - 2; pElemDOFNum[11] = _model.DOFMap[N4] * 3 - 1; pElemDOFNum[12] = _model.DOFMap[N4] * 3;
                    pElemDOFNum[13] = _model.DOFMap[N5] * 3 - 2; pElemDOFNum[14] = _model.DOFMap[N5] * 3 - 1; pElemDOFNum[15] = _model.DOFMap[N5] * 3;
                    pElemDOFNum[16] = _model.DOFMap[N6] * 3 - 2; pElemDOFNum[17] = _model.DOFMap[N6] * 3 - 1; pElemDOFNum[18] = _model.DOFMap[N6] * 3;
                    pElemDOFNum[19] = _model.DOFMap[N7] * 3 - 2; pElemDOFNum[20] = _model.DOFMap[N7] * 3 - 1; pElemDOFNum[21] = _model.DOFMap[N7] * 3;
                    pElemDOFNum[22] = _model.DOFMap[N8] * 3 - 2; pElemDOFNum[23] = _model.DOFMap[N8] * 3 - 1; pElemDOFNum[24] = _model.DOFMap[N8] * 3;
                    for i = 1:24
                        for j in bound
                            _RHS[pElemDOFNum[i]] -= (_K[i,j,_model.elemMatMap[e]]) * delta
                        end
                    end
                end
            end
        end
    elseif _model.rhsType == 0  # Domain
        println("   .Computing RHS - Domain!")
        for e = 1:_model.nElems
            N1 = 2 + (e - 1) % (nElemS) + div((e - 1) % (nElemS), _model.ny) + div(e - 1, nElemS) * nNodeS;
            N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
            N5 = N1 + nNodeS; N6 = N2 + nNodeS; N7 = N3 + nNodeS; N8 = N4 + nNodeS;
            pElemDOFNum[1]  = _model.DOFMap[N1] * 3 - 2; pElemDOFNum[2]  = _model.DOFMap[N1] * 3 - 1; pElemDOFNum[3]  = _model.DOFMap[N1] * 3;
            pElemDOFNum[4]  = _model.DOFMap[N2] * 3 - 2; pElemDOFNum[5]  = _model.DOFMap[N2] * 3 - 1; pElemDOFNum[6]  = _model.DOFMap[N2] * 3;
            pElemDOFNum[7]  = _model.DOFMap[N3] * 3 - 2; pElemDOFNum[8]  = _model.DOFMap[N3] * 3 - 1; pElemDOFNum[9]  = _model.DOFMap[N3] * 3;
            pElemDOFNum[10] = _model.DOFMap[N4] * 3 - 2; pElemDOFNum[11] = _model.DOFMap[N4] * 3 - 1; pElemDOFNum[12] = _model.DOFMap[N4] * 3;
            pElemDOFNum[13] = _model.DOFMap[N5] * 3 - 2; pElemDOFNum[14] = _model.DOFMap[N5] * 3 - 1; pElemDOFNum[15] = _model.DOFMap[N5] * 3;
            pElemDOFNum[16] = _model.DOFMap[N6] * 3 - 2; pElemDOFNum[17] = _model.DOFMap[N6] * 3 - 1; pElemDOFNum[18] = _model.DOFMap[N6] * 3;
            pElemDOFNum[19] = _model.DOFMap[N7] * 3 - 2; pElemDOFNum[20] = _model.DOFMap[N7] * 3 - 1; pElemDOFNum[21] = _model.DOFMap[N7] * 3;
            pElemDOFNum[22] = _model.DOFMap[N8] * 3 - 2; pElemDOFNum[23] = _model.DOFMap[N8] * 3 - 1; pElemDOFNum[24] = _model.DOFMap[N8] * 3;
            for i = 1:24
                _RHS[pElemDOFNum[i]] += _B[_axis + 1,i,_model.elemMatMap[e]]
            end
        end
    end
end

# Direct Solver: [K] 64 bits * nDOFs * nDOFs
function directMethod!(_model::Model, _x1::Vector{Float64}, _x2::Vector{Float64}, _x3::Vector{Float64}, _x4::Vector{Float64}, _x5::Vector{Float64}, _x6::Vector{Float64}, _RHS1::Vector{Float64}, _RHS2::Vector{Float64}, _RHS3::Vector{Float64}, _RHS4::Vector{Float64}, _RHS5::Vector{Float64}, _RHS6::Vector{Float64}, _K::Array{Float64,3})
    println("   .Direct Solver!")
    # Initializations
    K = spzeros(_model.nDOFs, _model.nDOFs)
    pElemDOFNum = zeros(UInt64, 24)
    N1::UInt64 = 0; N2::UInt64 = 0; N3::UInt64 = 0; N4::UInt64 = 0;
    N5::UInt64 = 0; N6::UInt64 = 0; N7::UInt64 = 0; N8::UInt64 = 0;  
    thisElemMat::UInt16 = 0  
    nElemS::UInt64 = _model.nx * _model.ny
    nNodeS::UInt64 = (_model.nx + 1) * (_model.ny + 1)
    # Assembly system matrix and solve for three rhs
    for e = 1:_model.nElems
        thisElemMat = _model.elemMatMap[e]
        N1 = 2 + (e - 1) % (nElemS) + div((e - 1) % (nElemS), _model.ny) + div(e - 1, nElemS) * nNodeS;
        N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
        N5 = N1 + nNodeS; N6 = N2 + nNodeS; N7 = N3 + nNodeS; N8 = N4 + nNodeS;
        pElemDOFNum[1]  = _model.DOFMap[N1] * 3 - 2; pElemDOFNum[2]  = _model.DOFMap[N1] * 3 - 1; pElemDOFNum[3]  = _model.DOFMap[N1] * 3;
        pElemDOFNum[4]  = _model.DOFMap[N2] * 3 - 2; pElemDOFNum[5]  = _model.DOFMap[N2] * 3 - 1; pElemDOFNum[6]  = _model.DOFMap[N2] * 3;
        pElemDOFNum[7]  = _model.DOFMap[N3] * 3 - 2; pElemDOFNum[8]  = _model.DOFMap[N3] * 3 - 1; pElemDOFNum[9]  = _model.DOFMap[N3] * 3;
        pElemDOFNum[10] = _model.DOFMap[N4] * 3 - 2; pElemDOFNum[11] = _model.DOFMap[N4] * 3 - 1; pElemDOFNum[12] = _model.DOFMap[N4] * 3;
        pElemDOFNum[13] = _model.DOFMap[N5] * 3 - 2; pElemDOFNum[14] = _model.DOFMap[N5] * 3 - 1; pElemDOFNum[15] = _model.DOFMap[N5] * 3;
        pElemDOFNum[16] = _model.DOFMap[N6] * 3 - 2; pElemDOFNum[17] = _model.DOFMap[N6] * 3 - 1; pElemDOFNum[18] = _model.DOFMap[N6] * 3;
        pElemDOFNum[19] = _model.DOFMap[N7] * 3 - 2; pElemDOFNum[20] = _model.DOFMap[N7] * 3 - 1; pElemDOFNum[21] = _model.DOFMap[N7] * 3;
        pElemDOFNum[22] = _model.DOFMap[N8] * 3 - 2; pElemDOFNum[23] = _model.DOFMap[N8] * 3 - 1; pElemDOFNum[24] = _model.DOFMap[N8] * 3;
        for i = 1:24
            for j = 1:24
                K[pElemDOFNum[i],pElemDOFNum[j]] += _K[i,j,thisElemMat];
            end
        end
    end
    # Solve for six RHS:
    _x1 .= K \ _RHS1
    _x2 .= K \ _RHS2
    _x3 .= K \ _RHS3
    _x4 .= K \ _RHS4
    _x5 .= K \ _RHS5
    _x6 .= K \ _RHS6
end

# Jacobi Preconditioner: assembly || M
function jacobiPrecond!(_model::Model, _M::Vector{Float64}, _K::Array{Float64,3})
    println("   .Jacobi Preconditioner!")
    # Initializations:
    pElemDOFNum = zeros(UInt64, 24)
    N1::UInt64 = 0; N2::UInt64 = 0; N3::UInt64 = 0; N4::UInt64 = 0;
    N5::UInt64 = 0; N6::UInt64 = 0; N7::UInt64 = 0; N8::UInt64 = 0;
    thisElemMat::UInt16 = 0
    nElemS::UInt64 = _model.nx * _model.ny
    nNodeS::UInt64 = (_model.nx + 1) * (_model.ny + 1)
    for e = 1:_model.nElems
        thisElemMat = _model.elemMatMap[e]
        N1 = 2 + (e - 1) % (nElemS) + div((e - 1) % (nElemS), _model.ny) + div(e - 1, nElemS) * nNodeS;
        N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
        N5 = N1 + nNodeS; N6 = N2 + nNodeS; N7 = N3 + nNodeS; N8 = N4 + nNodeS;
        pElemDOFNum[1]  = _model.DOFMap[N1] * 3 - 2; pElemDOFNum[2]  = _model.DOFMap[N1] * 3 - 1; pElemDOFNum[3]  = _model.DOFMap[N1] * 3;
        pElemDOFNum[4]  = _model.DOFMap[N2] * 3 - 2; pElemDOFNum[5]  = _model.DOFMap[N2] * 3 - 1; pElemDOFNum[6]  = _model.DOFMap[N2] * 3;
        pElemDOFNum[7]  = _model.DOFMap[N3] * 3 - 2; pElemDOFNum[8]  = _model.DOFMap[N3] * 3 - 1; pElemDOFNum[9]  = _model.DOFMap[N3] * 3;
        pElemDOFNum[10] = _model.DOFMap[N4] * 3 - 2; pElemDOFNum[11] = _model.DOFMap[N4] * 3 - 1; pElemDOFNum[12] = _model.DOFMap[N4] * 3;
        pElemDOFNum[13] = _model.DOFMap[N5] * 3 - 2; pElemDOFNum[14] = _model.DOFMap[N5] * 3 - 1; pElemDOFNum[15] = _model.DOFMap[N5] * 3;
        pElemDOFNum[16] = _model.DOFMap[N6] * 3 - 2; pElemDOFNum[17] = _model.DOFMap[N6] * 3 - 1; pElemDOFNum[18] = _model.DOFMap[N6] * 3;
        pElemDOFNum[19] = _model.DOFMap[N7] * 3 - 2; pElemDOFNum[20] = _model.DOFMap[N7] * 3 - 1; pElemDOFNum[21] = _model.DOFMap[N7] * 3;
        pElemDOFNum[22] = _model.DOFMap[N8] * 3 - 2; pElemDOFNum[23] = _model.DOFMap[N8] * 3 - 1; pElemDOFNum[24] = _model.DOFMap[N8] * 3;
        for i = 1:24
            _M[pElemDOFNum[i]] += _K[i,i,thisElemMat]
        end
    end
    _M .= _M .\ 1
end

# Preconditioned Conjugate Gradient Method:
function pcg_old!(_model::Model, _x::Vector{Float64}, _r::Vector{Float64}, _M::Vector{Float64}, _K::Array{Float64,3})
    println("   .PCG Solver!")
    # Initializations
    d = zeros(Float64, _model.nDOFs);
    q = zeros(Float64, _model.nDOFs);
    pElemDOFNum = zeros(UInt64, 24);
    N1 = 0; N2 = 0; N3 = 0; N4 = 0;
    N5 = 0; N6 = 0; N7 = 0; N8 = 0;
    q_temp = 0;
    nElemS = _model.nx * _model.ny; 
    nNodeS = (_model.nx + 1) * (_model.ny + 1);
    # PCG Initialization:
    d .= _r;
    d .*= _M;
    delta_new = (_r' * d)[1,1];
    delta_0 = delta_new;
    i_max = _model.pcgIter; 
    ii = 0;
    # PCG Iterations:
    while (ii < i_max) && (abs(delta_new) > _model.pcgTol * _model.pcgTol * abs(delta_0)) # (maximum(abs.(_r))>_pcgTol)
        @fastmath @inbounds @simd for e = 1:_model.nElems
            N1 = 2 + (e - 1) % (nElemS) + div((e - 1) % (nElemS), _model.ny) + div(e - 1, nElemS) * nNodeS;
            N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
            N5 = N1 + nNodeS; N6 = N2 + nNodeS; N7 = N3 + nNodeS; N8 = N4 + nNodeS;
            pElemDOFNum[1]  = _model.DOFMap[N1] * 3 - 2; pElemDOFNum[2]  = _model.DOFMap[N1] * 3 - 1; pElemDOFNum[3]  = _model.DOFMap[N1] * 3;
            pElemDOFNum[4]  = _model.DOFMap[N2] * 3 - 2; pElemDOFNum[5]  = _model.DOFMap[N2] * 3 - 1; pElemDOFNum[6]  = _model.DOFMap[N2] * 3;
            pElemDOFNum[7]  = _model.DOFMap[N3] * 3 - 2; pElemDOFNum[8]  = _model.DOFMap[N3] * 3 - 1; pElemDOFNum[9]  = _model.DOFMap[N3] * 3;
            pElemDOFNum[10] = _model.DOFMap[N4] * 3 - 2; pElemDOFNum[11] = _model.DOFMap[N4] * 3 - 1; pElemDOFNum[12] = _model.DOFMap[N4] * 3;
            pElemDOFNum[13] = _model.DOFMap[N5] * 3 - 2; pElemDOFNum[14] = _model.DOFMap[N5] * 3 - 1; pElemDOFNum[15] = _model.DOFMap[N5] * 3;
            pElemDOFNum[16] = _model.DOFMap[N6] * 3 - 2; pElemDOFNum[17] = _model.DOFMap[N6] * 3 - 1; pElemDOFNum[18] = _model.DOFMap[N6] * 3;
            pElemDOFNum[19] = _model.DOFMap[N7] * 3 - 2; pElemDOFNum[20] = _model.DOFMap[N7] * 3 - 1; pElemDOFNum[21] = _model.DOFMap[N7] * 3;
            pElemDOFNum[22] = _model.DOFMap[N8] * 3 - 2; pElemDOFNum[23] = _model.DOFMap[N8] * 3 - 1; pElemDOFNum[24] = _model.DOFMap[N8] * 3;
            for i = 1:24
                q_temp = 0;
                for j = 1:24
                    q_temp += _K[i,j,_model.elemMatMap[e]] * d[pElemDOFNum[j]];
                end
                q[pElemDOFNum[i]] += q_temp;
            end
        end
        alfa = delta_new / (d' * q)[1,1];
        d .*= alfa;
        _x .+= d;
        q .*= alfa;
        _r .-= q;
        q .= _r;
        q .*= _M;
        delta_old = delta_new;
        delta_new = (_r' * q)[1,1];
        beta = delta_new / delta_old;
        d .*= beta / alfa;
        d .+= q;
        q .*= 0;
        ii += 1;
    end
    println("    $ii steps")
    println("    Residue = ", sqrt(abs(delta_new) / abs(delta_0)))
end

# Preconditioned Conjugate Gradient Method:
function pcg!(_model::Model, _x::Vector{Float64}, _r::Vector{Float64}, _M::Vector{Float64}, _K::Array{Float64,3})
    println("   .PCG Solver!")
    # Initializations:    
    d = zeros(Float64, _model.nDOFs)
    q = zeros(Float64, _model.nDOFs)
    pElemDOFNum = zeros(UInt64,24)
    pElemDOFVar = zeros(Float64,24)
    thisElemMat::UInt16 = 0
    q_temp::Float64 = 0.0
    alfa::Float64 = 0.0
    beta::Float64 = 0.0
    ri::Float64   = 0.0
    qi::Float64   = 0.0
    delta_new::Float64 = 0.0
    delta_old::Float64 = 0.0
    N1::UInt64 = 0; N2::UInt64 = 0; N3::UInt64 = 0; N4::UInt64 = 0; N5::UInt64 = 0; N6::UInt64 = 0; N7::UInt64 = 0; N8::UInt64 = 0;
    nElem_perSlice::UInt64 = _model.nx * _model.ny
    nNode_perSlice::UInt64 = (_model.nx + 1) * (_model.ny + 1)
    # PCG Initialization:
    @inbounds for i=1:_model.nDOFs; d[i] = _r[i]*_M[i]; end
    delta_new = dot(_r,d)
    delta_0 = delta_new
    if (abs(delta_0)<1e-14); println("    x0 satisfied absolute tolerance criteria: delta < 1e-14"); return; end
    tolerance::Float64 = _model.pcgTol * _model.pcgTol * abs(delta_0);
    iteration_count::UInt64 = _model.pcgIter;
    # PCG Iterations:
    for ii = 1:_model.pcgIter
        # (EbE) q = Kd
        @fastmath @inbounds @simd for e = 1:_model.nElems
            thisElemMat = _model.elemMatMap[e];
            N1 = 2 + (e - 1) % (nElem_perSlice) + div((e - 1) % (nElem_perSlice), _model.ny) + div(e - 1, nElem_perSlice) * nNode_perSlice;
            N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
            N5 = N1 + nNode_perSlice; N6 = N2 + nNode_perSlice; N7 = N3 + nNode_perSlice; N8 = N4 + nNode_perSlice;
            pElemDOFNum[1]  = 3*_model.DOFMap[N1]-2; pElemDOFNum[2]  = pElemDOFNum[1] +1; pElemDOFNum[3]  = pElemDOFNum[2] +1; 
            pElemDOFNum[4]  = 3*_model.DOFMap[N2]-2; pElemDOFNum[5]  = pElemDOFNum[4] +1; pElemDOFNum[6]  = pElemDOFNum[5] +1; 
            pElemDOFNum[7]  = 3*_model.DOFMap[N3]-2; pElemDOFNum[8]  = pElemDOFNum[7] +1; pElemDOFNum[9]  = pElemDOFNum[8] +1; 
            pElemDOFNum[10] = 3*_model.DOFMap[N4]-2; pElemDOFNum[11] = pElemDOFNum[10]+1; pElemDOFNum[12] = pElemDOFNum[11]+1;
            pElemDOFNum[13] = 3*_model.DOFMap[N5]-2; pElemDOFNum[14] = pElemDOFNum[13]+1; pElemDOFNum[15] = pElemDOFNum[14]+1;
            pElemDOFNum[16] = 3*_model.DOFMap[N6]-2; pElemDOFNum[17] = pElemDOFNum[16]+1; pElemDOFNum[18] = pElemDOFNum[17]+1; 
            pElemDOFNum[19] = 3*_model.DOFMap[N7]-2; pElemDOFNum[20] = pElemDOFNum[19]+1; pElemDOFNum[21] = pElemDOFNum[20]+1; 
            pElemDOFNum[22] = 3*_model.DOFMap[N8]-2; pElemDOFNum[23] = pElemDOFNum[22]+1; pElemDOFNum[24] = pElemDOFNum[23]+1;
            @inbounds for i=1:24; pElemDOFVar[i] = d[pElemDOFNum[i]]; end
            @inbounds for i = 1:24
                q_temp = _K[i,1,thisElemMat] * pElemDOFVar[1]
                @inbounds for j=2:24; q_temp+= _K[i,j,thisElemMat] * pElemDOFVar[j]; end
                q[pElemDOFNum[i]] += q_temp
            end
        end
        alfa = delta_new / dot(d,q)
        delta_old = delta_new
        delta_new = 0.0
        @inbounds for i=1:_model.nDOFs
            _x[i] += d[i]*alfa
            _r[i] -= q[i]*alfa
            ri = _r[i]
            qi = ri*_M[i]
            q[i] = qi
            delta_new += ri*qi    
        end
        if (abs(delta_new) <= tolerance); iteration_count = ii; break; end
        beta = delta_new / delta_old
        @inbounds for i=1:_model.nDOFs
            d[i] *= beta
            d[i] += q[i]
            q[i] = 0.0;
        end
    end
    println("    $iteration_count steps");
    println("    Residue = ", sqrt(abs(delta_new) / abs(delta_0)));
end

# Compute Stress-FEM Effective property
function femEffective(_model::Model, _T::Vector{Float64}, _axis::Int, _B::Array{Float64,3})
    println("   .Updating Constitutive Matrix!")
    # Initialization
    SX::Float64 = 0.0; SY::Float64 = 0.0; SZ::Float64 = 0.0; SXY::Float64 = 0.0; SXZ::Float64 = 0.0; SYZ::Float64 = 0.0;
    N1::UInt64 = 0; N2::UInt64 = 0; N3::UInt64 = 0; N4::UInt64 = 0;
    N5::UInt64 = 0; N6::UInt64 = 0; N7::UInt64 = 0; N8::UInt64 = 0;
    nElemS::UInt64 = _model.nx * _model.ny
    nNodeS::UInt64 = (_model.nx + 1) * (_model.ny + 1)
    pElemDOFNum = zeros(UInt64, 24)
    bound = zeros(UInt8, 4)
    C = zeros(Float64, 6, 6)
    # Compute the effective properties for each test
    if _model.rhsType == 1    # Boundary
        deltaT::Float64 = 0.0
        if _axis == 0
            deltaT = _model.nx
            aux = 0
            bound[1] = 4; bound[2] = 7; bound[3] = 16; bound[4] = 19;
            for i = 1:_model.nz
                for eb = (nElemS + 1 - _model.ny + aux):(nElemS + aux)
                    for j in bound
                        SX  += (_B[1,j,_model.elemMatMap[eb]] * deltaT) 
                        SY  += (_B[2,j,_model.elemMatMap[eb]] * deltaT) 
                        SZ  += (_B[3,j,_model.elemMatMap[eb]] * deltaT)
                        SXY += (_B[4,j,_model.elemMatMap[eb]] * deltaT) 
                        SXZ += (_B[5,j,_model.elemMatMap[eb]] * deltaT) 
                        SYZ += (_B[6,j,_model.elemMatMap[eb]] * deltaT)
                    end
                end
                aux += nElemS
            end
        elseif _axis == 1
            deltaT = _model.ny
            bound[1] = 8; bound[2] = 11; bound[3] = 20; bound[4] = 23;
            for eb = 1:(_model.ny):_model.nElems
                for j in bound
                    SX  += (_B[1,j,_model.elemMatMap[eb]] * deltaT)
                    SY  += (_B[2,j,_model.elemMatMap[eb]] * deltaT) 
                    SZ  += (_B[3,j,_model.elemMatMap[eb]] * deltaT)
                    SXY += (_B[4,j,_model.elemMatMap[eb]] * deltaT) 
                    SXZ += (_B[5,j,_model.elemMatMap[eb]] * deltaT) 
                    SYZ += (_B[6,j,_model.elemMatMap[eb]] * deltaT)
                end
            end
        elseif _axis == 2
            deltaT = _model.nz
            bound[1] = 15; bound[2] = 18; bound[3] = 21; bound[4] = 24;
            for eb = (_model.nElems - nElemS + 1):_model.nElems
                for j in bound
                    SX  += (_B[1,j,_model.elemMatMap[eb]] * deltaT) 
                    SY  += (_B[2,j,_model.elemMatMap[eb]] * deltaT)
                    SZ  += (_B[3,j,_model.elemMatMap[eb]] * deltaT)
                    SXY += (_B[4,j,_model.elemMatMap[eb]] * deltaT)
                    SXZ += (_B[5,j,_model.elemMatMap[eb]] * deltaT)
                    SYZ += (_B[6,j,_model.elemMatMap[eb]] * deltaT)
                end
            end
        elseif _axis == 3
            deltaT = _model.ny
            bound[1] = 7; bound[2] = 10; bound[3] = 19; bound[4] = 22;
            for eb = 1:(_model.ny):_model.nElems
                for j in bound
                    SX  += (_B[1,j,_model.elemMatMap[eb]] * deltaT)
                    SY  += (_B[2,j,_model.elemMatMap[eb]] * deltaT) 
                    SZ  += (_B[3,j,_model.elemMatMap[eb]] * deltaT)
                    SXY += (_B[4,j,_model.elemMatMap[eb]] * deltaT)
                    SXZ += (_B[5,j,_model.elemMatMap[eb]] * deltaT)
                    SYZ += (_B[6,j,_model.elemMatMap[eb]] * deltaT)
                end
            end
        elseif _axis == 4
            deltaT = _model.nz
            bound[1] = 13; bound[2] = 16; bound[3] = 19; bound[4] = 22;
            for eb = (_model.nElems - nElemS + 1):_model.nElems
                for j in bound
                    SX  += (_B[1,j,_model.elemMatMap[eb]] * deltaT) 
                    SY  += (_B[2,j,_model.elemMatMap[eb]] * deltaT) 
                    SZ  += (_B[3,j,_model.elemMatMap[eb]] * deltaT)
                    SXY += (_B[4,j,_model.elemMatMap[eb]] * deltaT)
                    SXZ += (_B[5,j,_model.elemMatMap[eb]] * deltaT)
                    SYZ += (_B[6,j,_model.elemMatMap[eb]] * deltaT)
                end
            end
        elseif _axis == 5
            deltaT = _model.ny
            bound[1] = 9; bound[2] = 12; bound[3] = 21; bound[4] = 24;
            for eb = 1:(_model.ny):_model.nElems
                for j in bound
                    SX  += (_B[1,j,_model.elemMatMap[eb]] * deltaT)
                    SY  += (_B[2,j,_model.elemMatMap[eb]] * deltaT)
                    SZ  += (_B[3,j,_model.elemMatMap[eb]] * deltaT)
                    SXY += (_B[4,j,_model.elemMatMap[eb]] * deltaT)
                    SXZ += (_B[5,j,_model.elemMatMap[eb]] * deltaT)
                    SYZ += (_B[6,j,_model.elemMatMap[eb]] * deltaT)
                end
            end
        end
        for e = 1:_model.nElems
            N1 = 2 + (e - 1) % (nElemS) + div((e - 1) % (nElemS), _model.ny) + div(e - 1, nElemS) * nNodeS;
            N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
            N5 = N1 + nNodeS; N6 = N2 + nNodeS; N7 = N3 + nNodeS; N8 = N4 + nNodeS;
            pElemDOFNum[1]  = _model.DOFMap[N1] * 3 - 2; pElemDOFNum[2]  = _model.DOFMap[N1] * 3 - 1; pElemDOFNum[3]  = _model.DOFMap[N1] * 3;
            pElemDOFNum[4]  = _model.DOFMap[N2] * 3 - 2; pElemDOFNum[5]  = _model.DOFMap[N2] * 3 - 1; pElemDOFNum[6]  = _model.DOFMap[N2] * 3;
            pElemDOFNum[7]  = _model.DOFMap[N3] * 3 - 2; pElemDOFNum[8]  = _model.DOFMap[N3] * 3 - 1; pElemDOFNum[9]  = _model.DOFMap[N3] * 3;
            pElemDOFNum[10] = _model.DOFMap[N4] * 3 - 2; pElemDOFNum[11] = _model.DOFMap[N4] * 3 - 1; pElemDOFNum[12] = _model.DOFMap[N4] * 3;
            pElemDOFNum[13] = _model.DOFMap[N5] * 3 - 2; pElemDOFNum[14] = _model.DOFMap[N5] * 3 - 1; pElemDOFNum[15] = _model.DOFMap[N5] * 3;
            pElemDOFNum[16] = _model.DOFMap[N6] * 3 - 2; pElemDOFNum[17] = _model.DOFMap[N6] * 3 - 1; pElemDOFNum[18] = _model.DOFMap[N6] * 3;
            pElemDOFNum[19] = _model.DOFMap[N7] * 3 - 2; pElemDOFNum[20] = _model.DOFMap[N7] * 3 - 1; pElemDOFNum[21] = _model.DOFMap[N7] * 3;
            pElemDOFNum[22] = _model.DOFMap[N8] * 3 - 2; pElemDOFNum[23] = _model.DOFMap[N8] * 3 - 1; pElemDOFNum[24] = _model.DOFMap[N8] * 3;
            for i = 1:24
                SX  += (_B[1,i,_model.elemMatMap[e]] * _T[pElemDOFNum[i]])
                SY  += (_B[2,i,_model.elemMatMap[e]] * _T[pElemDOFNum[i]])
                SZ  += (_B[3,i,_model.elemMatMap[e]] * _T[pElemDOFNum[i]])
                SXY += (_B[4,i,_model.elemMatMap[e]] * _T[pElemDOFNum[i]])
                SXZ += (_B[5,i,_model.elemMatMap[e]] * _T[pElemDOFNum[i]])
                SYZ += (_B[6,i,_model.elemMatMap[e]] * _T[pElemDOFNum[i]])
            end
        end
    elseif _model.rhsType == 0 # Domain
        t = zeros(Float64, 24)
        if     (_axis == 0); t[4] = 1;  t[7] = 1;  t[16] = 1; t[19] = 1;
        elseif (_axis == 1); t[8] = 1;  t[11] = 1; t[20] = 1; t[23] = 1; 
        elseif (_axis == 2); t[15] = 1; t[18] = 1; t[21] = 1; t[24] = 1; 
        elseif (_axis == 3); t[7] = 1;  t[10] = 1; t[19] = 1; t[22] = 1; 
        elseif (_axis == 4); t[13] = 1; t[16] = 1; t[19] = 1; t[22] = 1; 
        elseif (_axis == 5); t[9] = 1;  t[12] = 1; t[21] = 1; t[24] = 1; end  
        for e = 1:_model.nElems
            N1 = 2 + (e - 1) % (nElemS) + div((e - 1) % (nElemS), _model.ny) + div(e - 1, nElemS) * nNodeS;
            N3 = N1 + _model.ny; N2 = N3 + 1; N4 = N1 - 1;
            N5 = N1 + nNodeS; N6 = N2 + nNodeS; N7 = N3 + nNodeS; N8 = N4 + nNodeS;
            pElemDOFNum[1]  = _model.DOFMap[N1] * 3 - 2; pElemDOFNum[2]  = _model.DOFMap[N1] * 3 - 1; pElemDOFNum[3]  = _model.DOFMap[N1] * 3;
            pElemDOFNum[4]  = _model.DOFMap[N2] * 3 - 2; pElemDOFNum[5]  = _model.DOFMap[N2] * 3 - 1; pElemDOFNum[6]  = _model.DOFMap[N2] * 3;
            pElemDOFNum[7]  = _model.DOFMap[N3] * 3 - 2; pElemDOFNum[8]  = _model.DOFMap[N3] * 3 - 1; pElemDOFNum[9]  = _model.DOFMap[N3] * 3;
            pElemDOFNum[10] = _model.DOFMap[N4] * 3 - 2; pElemDOFNum[11] = _model.DOFMap[N4] * 3 - 1; pElemDOFNum[12] = _model.DOFMap[N4] * 3;
            pElemDOFNum[13] = _model.DOFMap[N5] * 3 - 2; pElemDOFNum[14] = _model.DOFMap[N5] * 3 - 1; pElemDOFNum[15] = _model.DOFMap[N5] * 3;
            pElemDOFNum[16] = _model.DOFMap[N6] * 3 - 2; pElemDOFNum[17] = _model.DOFMap[N6] * 3 - 1; pElemDOFNum[18] = _model.DOFMap[N6] * 3;
            pElemDOFNum[19] = _model.DOFMap[N7] * 3 - 2; pElemDOFNum[20] = _model.DOFMap[N7] * 3 - 1; pElemDOFNum[21] = _model.DOFMap[N7] * 3;
            pElemDOFNum[22] = _model.DOFMap[N8] * 3 - 2; pElemDOFNum[23] = _model.DOFMap[N8] * 3 - 1; pElemDOFNum[24] = _model.DOFMap[N8] * 3;           
            for i = 1:24
                SX  += (_B[1,i,_model.elemMatMap[e]] * (t[i] - _T[pElemDOFNum[i]]))
                SY  += (_B[2,i,_model.elemMatMap[e]] * (t[i] - _T[pElemDOFNum[i]]))
                SZ  += (_B[3,i,_model.elemMatMap[e]] * (t[i] - _T[pElemDOFNum[i]]))
                SXY += (_B[4,i,_model.elemMatMap[e]] * (t[i] - _T[pElemDOFNum[i]]))
                SXZ += (_B[5,i,_model.elemMatMap[e]] * (t[i] - _T[pElemDOFNum[i]]))
                SYZ += (_B[6,i,_model.elemMatMap[e]] * (t[i] - _T[pElemDOFNum[i]]))
            end
        end
    end    
    C[_axis + 1,1] = SX / _model.nElems;  C[_axis + 1,2] = SY / _model.nElems;  C[_axis + 1,3] = SZ / _model.nElems;
    C[_axis + 1,4] = SXY / _model.nElems; C[_axis + 1,5] = SXZ / _model.nElems; C[_axis + 1,6] = SYZ / _model.nElems;
    return C
end

# Compute Effective Property
function computeProp(_C::Array{Float64,2})
    println(".Compute Effective Property!")
    S = inv(_C[1:3,1:3])
    E1 = 1.0 / S[1,1]; E2 = 1.0 / S[2,2]; E3 = 1.0 / S[3,3];
    G12 = _C[4,4];   G23 = _C[5,5];   G31 = _C[6,6];
    v21 = -S[1,2] / S[2,2]; v31 = -S[1,3] / S[3,3]; v32 = -S[2,3] / S[3,3];
    v12 = -S[2,1] / S[1,1]; v13 = -S[3,1] / S[1,1]; v23 = -S[3,2] / S[2,2];    
    println("$E1 E1");   println("$E2 E2");   println("$E3 E3");
    println("$v12 v12"); println("$v13 v13"); println("$v23 v23");
    println("$v21 v21"); println("$v31 v31"); println("$v32 v32");
    println("$G12 G12"); println("$G23 G23"); println("$G31 G31");
end

# -----------------
function homogenize(_arg)
    println("---------------------------")
    # Build the Model data struct:
    m_model = buildModel(_arg * ".json", _arg * ".raw")
    # Estimate Memory Consumption:
    estimateMemory(m_model)
    # SOLVE:
    println(".Solving")
    # Compute the stiffness matrix for each Material:
    m_K = zeros(Float64, 24, 24, m_model.nMat)
    m_B = zeros(Float64, 6, 24, m_model.nMat)
    m_C = zeros(Float64, 6, 6)
    elementStiffnessMatrices!(m_model, m_K, m_B)
    if (m_model.solverType == 0) # Preconditioned Conjugate Gradient Method
        # Initialize the effective tensor, the right hand side, the inicial guess and the preconditioner:  
        m_RHS = zeros(Float64, m_model.nDOFs)  
        m_X = zeros(Float64, m_model.nDOFs)
        m_M = zeros(Float64, m_model.nDOFs)
        # Compute the Jacobi preconditioner: 
        jacobiPrecond!(m_model, m_M, m_K)
        for axis = 0:5  # axis 0 = X || axis 1 = Y || axis 2 = Z || axis 3 = XY || axis 4 = XZ || axis 5 = YZ
            println("\n  Case ", axis + 1)
            # Compute the RHS: Boundary or Domain, rhsType: Boundary = 0 || Domain = 1
            computeRHS!(m_model, m_RHS, axis, m_K, m_B)
            # Solver (to ensure optimal RAM usage we call GC before and after the PCGM):    
            GC.gc()
            #pcg_old!(m_model, m_X, m_RHS, m_M, m_K);
            pcg!(m_model, m_X, m_RHS, m_M, m_K)
            GC.gc()
            # Compute Effective Property:
            m_C .+= femEffective(m_model, m_X, axis, m_B)
            m_RHS .*= 0
            m_X .*= 0
        end
    elseif (m_model.solverType == 1) # Direct Method
        # Compute the RHS: Boundary or Domain, rhsType: Boundary = 0 || Domain = 1
        # axis 0 = X || axis 1 = Y || axis 2 = Z || axis 3 = XY || axis 4 = XZ || axis 5 = YZ
        m_RHS1 = zeros(Float64, m_model.nDOFs); m_RHS2 = zeros(Float64, m_model.nDOFs); m_RHS3 = zeros(Float64, m_model.nDOFs);
        m_RHS4 = zeros(Float64, m_model.nDOFs); m_RHS5 = zeros(Float64, m_model.nDOFs); m_RHS6 = zeros(Float64, m_model.nDOFs);
        computeRHS!(m_model, m_RHS1, 0, m_K, m_B)
        computeRHS!(m_model, m_RHS2, 1, m_K, m_B)
        computeRHS!(m_model, m_RHS3, 2, m_K, m_B)
        computeRHS!(m_model, m_RHS4, 3, m_K, m_B)
        computeRHS!(m_model, m_RHS5, 4, m_K, m_B)
        computeRHS!(m_model, m_RHS6, 5, m_K, m_B)
        # Solver
        m_X1 = zeros(Float64, m_model.nDOFs); m_X2 = zeros(Float64, m_model.nDOFs); m_X3 = zeros(Float64, m_model.nDOFs);
        m_X4 = zeros(Float64, m_model.nDOFs); m_X5 = zeros(Float64, m_model.nDOFs); m_X6 = zeros(Float64, m_model.nDOFs);
        directMethod!(m_model, m_X1, m_X2, m_X3, m_X4, m_X5, m_X6, m_RHS1, m_RHS2, m_RHS3, m_RHS4, m_RHS5, m_RHS6, m_K)
        m_RHS1 = nothing; m_RHS2 = nothing; m_RHS3 = nothing;
        m_RHS4 = nothing; m_RHS5 = nothing; m_RHS6 = nothing;
        # Compute Effective Property:
        m_C .+= femEffective(m_model, m_X1, 0, m_B)
        m_C .+= femEffective(m_model, m_X2, 1, m_B)
        m_C .+= femEffective(m_model, m_X3, 2, m_B)
        m_C .+= femEffective(m_model, m_X4, 3, m_B)
        m_C .+= femEffective(m_model, m_X5, 4, m_B)
        m_C .+= femEffective(m_model, m_X6, 5, m_B)
        m_X1 = nothing; m_X2 = nothing; m_X3 = nothing;
        m_X4 = nothing; m_X5 = nothing; m_X6 = nothing;
    end
    println("---------------------------")
    println("Effective Properties:\n")
    for i = 1:6
        println("C[$i,:] = ", m_C[i,:])
    end
    println("\n--------------------------------------")
    return m_C
end

homogenize("benchmark1_100")
#homogenize("GGG40_300")

