def append_anti_basis_error(circuit, targets, p, basis):
    if p > 0:
        if basis == "X":
            circuit.append("Z_ERROR", targets, p)
        else:
            circuit.append("X_ERROR", targets, p)


class CircuitGenParameters:
    def __init__(
        self,
        rounds,
        after_clifford_depolarization=0,
        before_round_data_depolarization=0,
        before_measure_flip_probability=0,
        after_reset_flip_probability=0
    ):
        self.rounds = rounds
        self.after_clifford_depolarization = after_clifford_depolarization
        self.before_round_data_depolarization = before_round_data_depolarization
        self.before_measure_flip_probability = before_measure_flip_probability
        self.after_reset_flip_probability = after_reset_flip_probability
        self.burst_errors_depolarization = 0.5

    def append_begin_round_tick(self, circuit, data_qubits, ano_qubits):
        circuit.append("TICK")
        if self.before_round_data_depolarization > 0:
            norm_targets = []
            ano_targets = []
            for q in data_qubits:
                if q in ano_qubits:
                    ano_targets.append(q)
                else:
                    norm_targets.append(q)
            if norm_targets:
                circuit.append("DEPOLARIZE1", norm_targets, self.before_round_data_depolarization)
            if ano_targets:
                circuit.append("DEPOLARIZE1", ano_targets, self.burst_errors_depolarization)

    def append_unitary_1(self, circuit, name, targets, ano_qubits):
        circuit.append(name, targets)
        if self.after_clifford_depolarization > 0:
            norm_targets = []
            ano_targets = []
            for q in targets:
                if q in ano_qubits:
                    ano_targets.append(q)
                else:
                    norm_targets.append(q)
            if norm_targets:
                circuit.append("DEPOLARIZE1", norm_targets, self.after_clifford_depolarization)
            if ano_targets:
                circuit.append("DEPOLARIZE1", ano_targets, self.burst_errors_depolarization)

    def append_unitary_2(self, circuit, name, targets, ano_qubits):
        circuit.append(name, targets)
        if self.after_clifford_depolarization > 0:
            norm_targets = []
            ano_targets = []
            for i in range(0, len(targets), 2):
                if targets[i] in ano_qubits or targets[i+1] in ano_qubits:
                    ano_targets.extend(targets[i:i+2])
                else:
                    norm_targets.extend(targets[i:i+2])
            if norm_targets:
                circuit.append("DEPOLARIZE2", norm_targets, self.after_clifford_depolarization)
            if ano_targets:
                circuit.append("DEPOLARIZE2", ano_targets, self.burst_errors_depolarization)

    def append_reset(self, circuit, targets, basis="Z"):
        gate = "R" + basis
        circuit.append(gate, targets)
        append_anti_basis_error(circuit, targets, self.after_reset_flip_probability, basis)

    def append_measure(self, circuit, targets, basis="Z"):
        gate = "M" + basis
        append_anti_basis_error(circuit, targets, self.before_measure_flip_probability, basis)
        circuit.append(gate, targets)

    def append_measure_reset(self, circuit, targets, basis="Z"):
        gate = "MR" + basis
        append_anti_basis_error(circuit, targets, self.before_measure_flip_probability, basis)
        circuit.append(gate, targets)
        append_anti_basis_error(circuit, targets, self.after_reset_flip_probability, basis)

