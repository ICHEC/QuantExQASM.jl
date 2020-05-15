using SymPy

function ry_mat(ang)
    return [SymPy.cos(ang/2) -SymPy.sin(ang/2); SymPy.sin(ang/2) SymPy.cos(ang/2)]
end

function rz_mat(ang)
    return [SymPy.exp(-1im*ang/2) 0; 0 SymPy.exp(1im*ang/2)]
end

function gen_mat(a,b,c,d)
    return SymPy.exp(1im*d).*ry_mat(c)*rz_mat(b)*ry_mat(a)
end

function _params_zyz(mat::Matrix{<:Symbol})
    """Return the euler angles and phase for the ZYZ basis."""
    # We rescale the input matrix to be special unitary (det(U) = 1)
    # This ensures that the quaternion representation is real
    coeff = sympy.sqrt(sympy.det(mat))
    phase = -sympy.atan2(real(coeff),imag(coeff))
    su_mat = coeff * mat  # U in SU(2)
    # OpenQASM SU(2) parameterization:
    # U[0, 0] = exp(-i(phi+lambda)/2) * cos(theta/2)
    # U[0, 1] = -exp(-i(phi-lambda)/2) * sin(theta/2)
    # U[1, 0] = exp(i(phi-lambda)/2) * sin(theta/2)
    # U[1, 1] = exp(i(phi+lambda)/2) * cos(theta/2)
    theta = 2 * sympy.atan2(abs(su_mat[2, 1]), abs(su_mat[1, 1]))
    phiplambda = 2 * sympy.atan2(real(su_mat[2, 2]), imag(su_mat[2, 2]))
    phimlambda = 2 * sympy.atan2(real(su_mat[2, 1]), imag(su_mat[2, 1]))
    phi = (phiplambda + phimlambda) / 2
    lam = (phiplambda - phimlambda) / 2
    return theta, phi, lam, phase
end

function _params_zyz(mat::Matrix{<:Number})
    """Return the euler angles and phase for the ZYZ basis."""
    # We rescale the input matrix to be special unitary (det(U) = 1)
    # This ensures that the quaternion representation is real
    coeff = sqrt(det(mat) +0im)
    phase = -atan(real(coeff), imag(coeff))
    su_mat = coeff * mat  # U in SU(2)
    # OpenQASM SU(2) parameterization:
    # U[0, 0] = exp(-i(phi+lambda)/2) * cos(theta/2)
    # U[0, 1] = -exp(-i(phi-lambda)/2) * sin(theta/2)
    # U[1, 0] = exp(i(phi-lambda)/2) * sin(theta/2)
    # U[1, 1] = exp(i(phi+lambda)/2) * cos(theta/2)
    theta = 2 * atan(abs(su_mat[2, 1]), abs(su_mat[1, 1]))
    phiplambda = 2 * atan(real(su_mat[2, 2]), imag(su_mat[2, 2]))
    phimlambda = 2 * atan(real(su_mat[2, 1]), imag(su_mat[2, 1]))
    phi = (phiplambda + phimlambda) / 2
    lam = (phiplambda - phimlambda) / 2
    return theta, phi, lam, phase
end