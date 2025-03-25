# this is very slow! change this function to remove the Tensors and TensorOperations
# packages and write it by hand!
function get_zem_stiffness!(storage, params, Kinv, ::NoRotation, i)
    (; C_1) = storage
    (; C) = params
    C_1 .= 0.0

    # Compute C_1 from C and Kinv
    for i in 1:3, j in 1:3
        sum_value = 0.0
        for k in 1:3, l in 1:3
            sum_value += C[i, j, k, l] * Kinv[k, l]
        end
        C_1[i, j] = sum_value
    end

    return C_1
end

# this is very slow! change this function to remove the Tensors and TensorOperations
# packages and write it by hand!
function get_zem_stiffness!(storage, params, Kinv, ::FlanaganTaylorRotation, i)
    (; C_1, C_rotated) = storage
    (; C) = params
    C_1 .= 0.0
    C_rotated .= 0.0
    R = get_tensor(storage.rotation, i)

    # Rotate the C tensor
    for m in 1:3, n in 1:3, o in 1:3, p in 1:3
        sum_value = 0.0
        for i in 1:3, j in 1:3, k in 1:3, l in 1:3
            sum_value += C[i, j, k, l] * R[i, m] * R[j, n] * R[k, o] * R[l, p]
        end
        C_rotated[m, n, o, p] = sum_value
    end

    # Compute C_1 from C_rotated and Kinv
    for i in 1:3, j in 1:3
        sum_value = 0.0
        for k in 1:3, l in 1:3
            sum_value += C_rotated[i, j, k, l] * Kinv[k, l]
        end
        C_1[i, j] = sum_value
    end

    return C_1
end

function get_hooke_matrix(nu, λ, μ)
    a = (1 - nu) * λ / nu
    CVoigt = SMatrix{6,6,Float64,36}(
        a, λ, λ, 0, 0, 0,
        λ, a, λ, 0, 0, 0,
        λ, λ, a, 0, 0, 0,
        0, 0, 0, μ, 0, 0,
        0, 0, 0, 0, μ, 0,
        0, 0, 0, 0, 0, μ
    )

    C = MArray{Tuple{3,3,3,3}, Float64}(undef)

    # Voigt notation index mapping
    voigt_map = @SVector [(1,1), (2,2), (3,3), (2,3), (1,3), (1,2)]

    # Convert Voigt matrix to full 3x3x3x3 tensor
    for I in 1:6, J in 1:6
        i1, i2 = voigt_map[I]
        j1, j2 = voigt_map[J]
        C[i1, i2, j1, j2] = CVoigt[I, J]
        C[i2, i1, j1, j2] = CVoigt[I, J]
        C[i1, i2, j2, j1] = CVoigt[I, J]
        C[i2, i1, j2, j1] = CVoigt[I, J]
    end

    return C
end
