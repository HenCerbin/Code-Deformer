from circuit_gen_params import CircuitGenParameters
from code_deformation import LogicalQubit
from itertools import chain
import stim


class MeasurementRecord:
    def __init__(self):
        self.t = 0
        self.record = {}

    def measure(self, measurements):
        for i, measurement in enumerate(measurements):
            try:
                self.record[measurement].append(self.t + i)
            except KeyError:
                self.record[measurement] = [self.t + i]
        self.t += len(measurements)

    def measure_rec(self, measurement, idx):
        return stim.target_rec(self.record[measurement][idx] - self.t)


def _generate_unshell_surface_code_circuit(params, logical_qubit, is_memory_z: bool):
    if params.rounds < 1:
        raise AttributeError("Need rounds >= 1.")

    chosen_basis = "XZ"[is_memory_z]
    data_coords = logical_qubit.data_coords
    ano_coords = logical_qubit.ano_coords
    stab_coords = {basis: set(coord for coord in logical_qubit.stabs[basis].keys()) for basis in "XZ"}
    gauge_coords = {basis: set(coord for coord in logical_qubit.gauges[basis].keys()) for basis in "XZ"}

    # Index the measurement qubits and data qubits.
    p2q = {}
    for p in chain(logical_qubit.qubit_coords, logical_qubit.ano_coords):
        p2q[p] = logical_qubit.coord_to_index(p)

    # Make target lists for various types of qubits.
    data_qubits = sorted([p2q[p] for p in data_coords])
    ano_qubits = sorted([p2q[p] for p in ano_coords])
    stab_qubits = {basis: sorted([p2q[p] for p in stab_coords[basis]]) for basis in "XZ"}
    gauge_ancilla_qubits = {basis: sorted([p2q[p] for p in gauge_coords[basis] if p not in data_coords]) for basis in "XZ"}
    gauge_data_qubits = {basis: sorted([p2q[p] for p in gauge_coords[basis] if p in data_coords]) for basis in "XZ"}

    # List out CNOT gate targets using given interaction orders.
    order = {"X": [(1, 1), (-1, 1), (1, -1), (-1, -1)],
             "Z": [(1, 1), (1, -1), (-1, 1), (-1, -1)]}
    stab_cnot_targets = [[] for _ in range(4)]
    gauge_cnot_targets = {basis: [[] for _ in range(4)] for basis in "XZ"}
    for k in range(4):
        for basis in "XZ":
            for coord in logical_qubit.stabs[basis].keys():
                data = (coord[0] + order[basis][k][0], coord[1] + order[basis][k][1])
                if data in logical_qubit.stabs[basis][coord]:
                    stab_cnot_targets[k].append(p2q[coord if basis == "X" else data])
                    stab_cnot_targets[k].append(p2q[data if basis == "X" else coord])
            for coord in logical_qubit.gauges[basis].keys():
                data = (coord[0] + order[basis][k][0], coord[1] + order[basis][k][1])
                if data in logical_qubit.gauges[basis][coord]:
                    gauge_cnot_targets[basis][k].append(p2q[coord if basis == "X" else data])
                    gauge_cnot_targets[basis][k].append(p2q[data if basis == "X" else coord])

    # Build the repeated actions that make up the surface code cycle.
    record = MeasurementRecord()

    def _generate_cycle_actions(is_gauge_z):
        cycle_actions = stim.Circuit()
        gauge_basis = "XZ"[is_gauge_z]
        x_qubits = stab_qubits["X"] + ([] if is_gauge_z else gauge_ancilla_qubits["X"])
        cycle_actions.append("TICK")
        params.append_reset(cycle_actions, stab_qubits["X"] + stab_qubits["Z"] + gauge_ancilla_qubits[gauge_basis])
        params.append_begin_round_tick(cycle_actions, data_qubits, ano_qubits)
        params.append_unitary_1(cycle_actions, "H", x_qubits, ano_qubits)
        for k in range(4):
            cycle_actions.append("TICK")
            params.append_unitary_2(cycle_actions, "CNOT",
                                    stab_cnot_targets[k] + gauge_cnot_targets[gauge_basis][k], ano_qubits)
        cycle_actions.append("TICK")
        params.append_unitary_1(cycle_actions, "H",
                                x_qubits + ([] if is_gauge_z else gauge_data_qubits["X"]), ano_qubits)
        cycle_actions.append("TICK")
        params.append_measure(cycle_actions, stab_qubits["X"] + stab_qubits["Z"] +
                              gauge_ancilla_qubits[gauge_basis] + gauge_data_qubits[gauge_basis])
        if not is_gauge_z and gauge_data_qubits["X"]:
            params.append_unitary_1(cycle_actions, "H", gauge_data_qubits["X"], ano_qubits)

        record.measure([("X", "stab", qubit) for qubit in stab_qubits["X"]] +
                       [("Z", "stab", qubit) for qubit in stab_qubits["Z"]] +
                       [(gauge_basis, "gauge", qubit) for qubit in gauge_ancilla_qubits[gauge_basis]] +
                       [(gauge_basis, "gauge", qubit) for qubit in gauge_data_qubits[gauge_basis]])

        return cycle_actions

    def _generate_stab_detectors():
        detectors = stim.Circuit()
        for basis in "XZ":
            for coord in logical_qubit.stabs[basis].keys():
                detectors.append(
                    "DETECTOR",
                    [record.measure_rec((basis, "stab", p2q[coord]), -1),
                     record.measure_rec((basis, "stab", p2q[coord]), -2)],
                    coord + (0,)
                )
        return detectors

    def _generate_gauge_detectors(is_gauge_z):
        detectors = stim.Circuit()
        basis = "XZ"[is_gauge_z]
        for super_stab in logical_qubit.super_stabs[basis]:
            detectors.append(
                "DETECTOR",
                [record.measure_rec((basis, "gauge", p2q[coord]), -1) for coord in super_stab] +
                [record.measure_rec((basis, "gauge", p2q[coord]), -2) for coord in super_stab],
                (-1, -1, 0)
            )
        return detectors

    # Build the start of the circuit, getting a state that's ready to cycle.
    # In particular, the first cycle has different detectors and so has to be handled special.
    head = stim.Circuit()
    for p, q in sorted(p2q.items(), key=lambda item: item[1]):
        head.append("QUBIT_COORDS", q, p)
    params.append_reset(head, data_qubits, chosen_basis)
    head += _generate_cycle_actions(is_memory_z)
    for coord in logical_qubit.stabs[chosen_basis].keys():
        head.append(
            "DETECTOR",
            [record.measure_rec((chosen_basis, "stab", p2q[coord]), -1)],
            coord + (0,)
        )
    for super_stab in logical_qubit.super_stabs[chosen_basis]:
        head.append(
            "DETECTOR",
            [record.measure_rec((chosen_basis, "gauge", p2q[coord]), -1) for coord in super_stab],
            (-1, -1, 0)
        )
    head += _generate_cycle_actions(not is_memory_z)
    head.append("SHIFT_COORDS", [], (0, 0, 1))
    head += _generate_stab_detectors()

    # Build the repeated body of the circuit, including the detectors comparing to previous cycles.
    body = stim.Circuit()
    body += _generate_cycle_actions(is_memory_z)
    body.append("SHIFT_COORDS", [], (0, 0, 1))
    body += _generate_stab_detectors()
    body += _generate_gauge_detectors(is_memory_z)
    body += _generate_cycle_actions(not is_memory_z)
    body.append("SHIFT_COORDS", [], (0, 0, 1))
    body += _generate_stab_detectors()
    body += _generate_gauge_detectors(not is_memory_z)

    # Build the end of the circuit, getting out of the cycle state and terminating.
    # In particular, the data measurements create detectors that have to be handled special.
    # Also, the tail is responsible for identifying the logical observable.
    tail = stim.Circuit()
    params.append_measure(tail, data_qubits, chosen_basis)
    record.measure([(chosen_basis, "data", qubit) for qubit in data_qubits])
    # Detectors.
    for coord, acting_coords in logical_qubit.stabs[chosen_basis].items():
        tail.append(
            "DETECTOR",
            [record.measure_rec((chosen_basis, "data", p2q[act_coord]), -1) for act_coord in acting_coords] +
            [record.measure_rec((chosen_basis, "stab", p2q[coord]), -1)],
            coord + (1,)
        )
    for super_stab in logical_qubit.super_stabs[chosen_basis]:
        detector = []
        for coord in super_stab:
            for act_coord in logical_qubit.gauges[chosen_basis][coord]:
                detector.append(record.measure_rec((chosen_basis, "data", p2q[act_coord]), -1))
        detector += [record.measure_rec((chosen_basis, "gauge", p2q[coord]), -1) for coord in super_stab]
        tail.append("DETECTOR", detector, (-1, -1, 1))
    # Logical observable
    tail.append(
        "OBSERVABLE_INCLUDE",
        [record.measure_rec((chosen_basis, "data", p2q[coord]), -1) for coord in logical_qubit.observable[chosen_basis]],
        0)

    # Combine to form final circuit.
    full_circuit = head + body * (params.rounds - 1) + tail
    return full_circuit


def generate_surface_code_circuit(params, logical_qubit, is_shell, is_memory_z):
    if is_shell:
        pass
        # return _generate_shell_surface_code_circuit(params, logical_qubit, is_memory_x)
    else:
        return _generate_unshell_surface_code_circuit(params, logical_qubit, is_memory_z)


if __name__ == "__main__":
    d = 3
    is_shell = False
    is_memory_z = True
    noise = 0.001

    P = CircuitGenParameters(d, noise, noise, noise, noise)
    Q = LogicalQubit(d, True)

    # for i in range(6):
    #     for j in [4, 5]:
    #         Q.disable((i * 2 + 1, j * 2 + 1))
    #
    # Q.disable((6, 4))
    # Q.disable((1, 1))
    # Q.disable((1, 7))
    # Q.disable((5, 1))
    # # Q.disable((5, 7))
    # # Q.disable((7, 1))
    # # Q.disable((7, 7))
    # # Q.disable((11, 1))
    # # Q.disable((11, 7))
    Circuit = generate_surface_code_circuit(P, Q, False, True)
    pass
