module QFT

export gen_qft, gen_iqft
using ..GateOps

"""
    gen_qft(qubit_indices::Vector{Int})

Generates an OpenQASM representation of the quantum Fourier transform
for the given number of qubits.

"""
function gen_qft(qubit_indices::Vector{Int})
    circuit = Vector{GateOps.GateCall}()

    angle(k) = 2*pi / 2^k;
    for (iidx,ival) = Iterators.reverse(enumerate(qubit_indices))
        push!(circuit, GateOps.hadamard(qubit_indices[iidx]))
        for (jidx,jval) = Iterators.reverse(enumerate(qubit_indices[1:iidx-1]))
            theta = angle(iidx - jidx);
            push!(circuit, GateOps.c_r_phase(ival, jval, theta))
        end
    end
    for idx = 1:div(length(qubit_indices),2)
        push!(circuit, GateOps.swap(qubit_indices[idx], qubit_indices[length(qubit_indices)-idx+1]))
    end
    return circuit
end

"""
    gen_iqft(qubit_indices::Vector{Int})

Generates an OpenQASM representation of the inverse quantum Fourier transform
for the given number of qubits.

"""
function gen_iqft(qubit_indices::Vector{Int})
    circuit = Vector{GateOps.GateCall}()

    angle(k) = -2*pi / 2^k;
    for (iidx,ival) = enumerate(qubit_indices)
        for (jidx,jval) = enumerate(qubit_indices[1:iidx-1])
            theta = angle(iidx - jidx);
            push!(circuit, GateOps.c_r_phase(ival, jval, theta))
        end
        circuit *= GateOps.apply_gate_h(q_reg, qubit_indices[iidx])

    end
    for idx = 1:div(length(qubit_indices),2)
        push!(circuit, GateOps.swap(qubit_indices[idx], qubit_indices[length(qubit_indices)-idx+1]))
    end
    return circuit
end

end
