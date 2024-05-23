import networkx as nx
from itertools import chain, product


def anti_commute(stab1: set, stab2: set) -> bool:
    return len(stab1 & stab2) % 2 == 1


class LogicalQubit:
    def __init__(self, distance, is_rotated: bool):
        self.distance = {basis: distance for basis in "XZ"}
        self.qubit_coords = set()
        self.data_coords = set()
        self.ano_coords = set()
        self.defect_coords = set()

        self.stabs = {basis: {} for basis in "XZ"}
        self.gauges = {basis: {} for basis in "XZ"}
        self.super_stabs = {basis: [] for basis in "XZ"}

        self.observable = {basis: set() for basis in "XZ"}
        # Each qubit on edge[basis] has only one $basis$ stabilizer acting on it
        self.edges = {basis: [set(), set()] for basis in "XZ"}
        # corners[i][j] is the intersection of edge["X"][i] and edge["Z"][j]
        self.corners = [[None, None],
                        [None, None]]

        self.decode_graph = {basis: None for basis in "XZ"}

        if is_rotated:
            self.generate_rotated_surface_code()
            self.coord_to_index = lambda q: int(q[0] + (q[1] - q[0] % 2) * (distance + 0.5))
        else:
            pass

    def generate_rotated_surface_code(self):
        d = self.distance["X"]

        # Place data qubits.
        for x in range(d):
            for y in range(d):
                q = (x * 2 + 1, y * 2 + 1)
                self.data_coords.add(q)

                if x == 0:
                    self.observable["X"].add(q)
                if y == 0:
                    self.observable["Z"].add(q)

                if x == 0:
                    self.edges["X"][0].add(q)
                elif x == d - 1:
                    self.edges["X"][1].add(q)
                if y == 0:
                    self.edges["Z"][0].add(q)
                elif y == d - 1:
                    self.edges["Z"][1].add(q)

        self.corners = [[(1, 1), (1, d * 2 - 1)],
                        [(d * 2 - 1, 1), (d * 2 - 1, d * 2 - 1)]]

        # Place measurement qubits.
        def neighbors(q):
            qubits = set()
            for i in [-1, 1]:
                for j in [-1, 1]:
                    qubit = (q[0] + i, q[1] + j)
                    if qubit in self.data_coords:
                        qubits.add(qubit)
            return qubits

        for x in range(d + 1):
            for y in range(d + 1):
                q = (x * 2, y * 2)
                on_boundary_1 = x == 0 or x == d
                on_boundary_2 = y == 0 or y == d
                parity = x % 2 != y % 2
                if on_boundary_1 and parity:
                    pass
                elif on_boundary_2 and not parity:
                    pass
                elif parity:
                    self.stabs["X"][q] = neighbors(q)
                else:
                    self.stabs["Z"][q] = neighbors(q)

        self._check()

    def burst_error(self, coord, r):
        for q in self.qubit_coords:
            if pow(q[0] - coord[0], 2) + pow(q[1] - coord[1], 2) <= 2 * pow(r, 2):
                self.ano_coords.add(q)

    def super_stabilizer(self, basis, gauge_coords):
        super_stab = set()
        for coord in gauge_coords:
            super_stab.symmetric_difference_update(self.gauges[basis][coord])
        return super_stab

    def disable(self, coord):
        self.defect_coords.add(coord)
        if coord in self.data_coords:
            self._disable_data(coord)
        elif coord in self.qubit_coords:
            self._disable_ancilla(coord)

    def update_distance(self):
        if self.observable:
            for basis, basis2 in ["XZ", "ZX"]:
                G = nx.Graph()
                G.add_nodes_from(self.stabs[basis].keys())
                G.add_nodes_from(["s%d" % k for k in range(len(self.super_stabs[basis]))])
                G.add_nodes_from(["e0", "e1"])
                G_edges = {q: [] for q in self.data_coords}

                for coord, stab in self.stabs[basis].items():
                    for q in stab:
                        G_edges[q].append(coord)
                for idx, super_stab in enumerate(self.super_stabs[basis]):
                    for q in self.super_stabilizer(basis, super_stab):
                        G_edges[q].append("s%d" % idx)
                for k in range(2):
                    for q in self.edges[basis][k]:
                        G_edges[q].append("e%d" % k)

                for q in self.data_coords:
                    if G_edges[q]:
                        G.add_edge(*G_edges[q])
                self.distance[basis2] = nx.shortest_path_length(G, "e0", "e1")
        else:
            self.distance = {"X": 1, "Z": 1}

    def _disable_data(self, coord):
        # print(coord)
        while coord in self.data_coords:
            rel_edge_idx = {basis: -1 for basis in "XZ"}
            for basis, idx in product("XZ", range(2)):
                if coord in self.edges[basis][idx]:
                    rel_edge_idx[basis] = idx

            if rel_edge_idx["X"] == -1 and rel_edge_idx["Z"] == -1:
                for basis in "XZ":
                    self._add_gauge(basis, coord)
            else:
                anti_stab = {basis: set() for basis in "XZ"}
                for basis in "XZ":
                    for stab in self.stabs[basis].values():
                        if coord in stab:
                            anti_stab[basis] = stab
                    for gauge_coords in self.super_stabs[basis]:
                        stab = self.super_stabilizer(basis, gauge_coords)
                        if coord in stab:
                            anti_stab[basis] = stab

                if rel_edge_idx["X"] == -1 or rel_edge_idx["Z"] == -1:
                    # coord is on edge[basis]
                    # Measure coord on $basis2$ => Add $basis2$ gauge and fix it
                    # edge[basis] = edge[basis] + anti_stab[basis]
                    basis, basis2 = "XZ" if rel_edge_idx["X"] != -1 else "ZX"

                    # check anti_stab[basis]
                    if anti_stab[basis]:
                        for idx2 in range(2):
                            if anti_stab[basis].intersection(self.edges[basis2][idx2]):
                                rel_edge_idx[basis2] = idx2

                        if rel_edge_idx[basis2] != -1:
                            # if anti_stab[basis] also acts on edge[basis2], disable the corner
                            self._disable_data(self.corners[rel_edge_idx["X"]][rel_edge_idx["Z"]])
                        else:
                            self._add_gauge(basis2, coord)
                            self._fix_gauge(basis2, coord)
                    else:
                        # if unable to find anti_stab[basis], shrink edge[basis2] to make coord a corner
                        if len(self.edges[basis2][0]) < len(self.edges[basis2][1]):
                            rel_edge_idx[basis2] = 0
                        else:
                            rel_edge_idx[basis2] = 1
                        self._disable_data(self.corners[rel_edge_idx["X"]][rel_edge_idx["Z"]])

                else:
                    # coord is on corner of rel_edge["X"] and rel_edge["Z"]
                    rel_edge = {basis: self.edges[basis][rel_edge_idx[basis]] for basis in "XZ"}

                    # Determine $basis$ of the measurement
                    # edge[basis] -= {coord}
                    # edge[basis2] = edge[basis2] * anti_stab[basis2]
                    if not anti_stab["X"] or not anti_stab["Z"]:
                        basis = "X" if anti_stab["Z"] else "Z"
                    elif len(rel_edge["X"]) == len(rel_edge["Z"]):
                        basis = "X" if len(anti_stab["X"]) > len(anti_stab["Z"]) else "Z"
                    else:
                        basis = "X" if len(rel_edge["X"]) > len(rel_edge["Z"]) else "Z"
                    basis2 = "XZ"[basis == "X"]

                    def new_corner(q):
                        anti_stab_k = anti_stab[basis2].difference({q})

                        edge_segments = []
                        for measurement in chain(self.stabs[basis2].values(), self.gauges[basis2].values()):
                            edge_segment = measurement.intersection(edge)
                            if edge_segment:
                                edge_segments.append(edge_segment)

                        while edge.intersection(anti_stab_k):
                            loop_flag = False
                            for edge_segment in edge_segments:
                                if q in edge_segment:
                                    edge_segment.remove(q)
                                    q = edge_segment.pop()
                                    loop_flag = True
                                    break
                            anti_stab_k.discard(q)
                            assert loop_flag

                        return q

                    # check if anti_stab[basis2] connect both edge[basis]
                    edge = self.edges[basis][1 - rel_edge_idx[basis]]
                    i = 1 - rel_edge_idx["X"] if basis == "X" else rel_edge_idx["X"]
                    j = 1 - rel_edge_idx["Z"] if basis == "Z" else rel_edge_idx["Z"]
                    if anti_stab[basis2].intersection(edge):
                        if self.corners[i][j] not in anti_stab[basis2]:
                            # similar to the case if coord is on the edge[basis]
                            self._disable_data(self.corners[i][j])
                            continue
                        else:
                            self.corners[i][j] = new_corner(self.corners[i][j])

                    edge = self.edges[basis][rel_edge_idx[basis]]
                    self.corners[rel_edge_idx["X"]][rel_edge_idx["Z"]] = new_corner(coord)

                    self._add_gauge(basis, coord)
                    self._fix_gauge(basis, coord)

            self._check()

    def _disable_ancilla(self, coord):
        basis, basis2 = "XZ" if coord in chain(self.stabs["X"].keys(), self.gauges["X"].keys()) else "ZX"
        # If $basis$ stabilizer/gauge contains qubit in $basis2$ edge, the count is exactly 2
        while coord in chain(self.stabs[basis].keys(), self.gauges[basis].keys()):
            measurement = self.stabs[basis][coord] if coord in self.stabs[basis].keys() else self.gauges[basis][coord]
            super_stab = set()
            edge_qubits = set()
            corner_qubit = set()
            for q in measurement:
                if any(q in edge for edge in self.edges[basis2]):
                    edge_qubits.add(q)
                    if any(q in edge for edge in self.edges[basis]):
                        corner_qubit.add(q)

            if edge_qubits:
                # if len(edge_qubits) > 2:
                #     raise Exception("Too small distance")
                sum_coord = tuple(map(sum, zip(*edge_qubits)))
                coord2 = (sum_coord[0] - coord[0], sum_coord[1] - coord[1])
                if len(edge_qubits) == 2 and coord2 not in self.qubit_coords.union(self.defect_coords):
                    self.gauges[basis][coord2] = edge_qubits
                    super_stab.add(coord2)
                else:
                    self._disable_data(corner_qubit.pop() if corner_qubit else edge_qubits.pop())
                    continue

            for q in measurement.difference(edge_qubits):
                self._add_gauge(basis, q)
                super_stab.add(q)

            if coord in self.stabs[basis].keys():
                self.stabs[basis].pop(coord)
                self.super_stabs[basis].append(super_stab)
            elif coord in self.gauges[basis].keys():
                self.gauges[basis].pop(coord)
                for super_stab2 in self.super_stabs[basis]:
                    if coord in super_stab2:
                        super_stab2.remove(coord)
                        super_stab2.symmetric_difference_update(super_stab)

            self._check()

    def _add_gauge(self, basis, coord):
        if coord not in self.gauges[basis].keys():
            basis2 = "XZ"[basis == "X"]

            anti_stabs = []  # coordinates of anti-commute stabilizers
            for coord2, stab2 in self.stabs[basis2].items():
                if coord in stab2:
                    anti_stabs.append(coord2)

            anti_super_stabs = []  # indices of anti-commute super-stabilizers
            for idx, gauge_coords in enumerate(self.super_stabs[basis2]):
                if coord in self.super_stabilizer(basis2, gauge_coords):
                    anti_super_stabs.append(idx)

            # Information Preservation
            if anti_stabs:
                stab = self.stabs[basis2][anti_stabs[0]]
            elif anti_super_stabs:
                stab = self.super_stabilizer(basis2, self.super_stabs[basis2][anti_super_stabs[0]])
            else:
                stab = set()

            if coord in self.observable[basis2]:
                assert stab
                self.observable[basis2].symmetric_difference_update(stab)

            for edge in self.edges[basis2]:
                if coord in edge:
                    assert stab
                    edge.symmetric_difference_update(stab)

            # Add new gauges
            self.gauges[basis][coord] = {coord}
            for coord2 in anti_stabs:
                self.gauges[basis2][coord2] = self.stabs[basis2].pop(coord2)

            # Update super-stabilizers
            if len(anti_stabs) == 2:
                self.super_stabs[basis2].append(set(anti_stabs))
            elif len(anti_super_stabs) == 2:
                idx1, idx2 = anti_super_stabs
                self.super_stabs[basis2][idx1].symmetric_difference_update(self.super_stabs[basis2].pop(idx2))
            elif len(anti_stabs) == len(anti_super_stabs) == 1:
                self.super_stabs[basis2][anti_super_stabs[0]].add(anti_stabs[0])
            elif len(anti_super_stabs) == 1:
                self.super_stabs[basis2].pop(anti_super_stabs[0])

    def _fix_gauge(self, basis, coord):
        # Fix gauge
        gauge = self.gauges[basis].pop(coord)
        self.stabs[basis][coord] = gauge
        # Update super-stabilizers
        for super_stab in self.super_stabs[basis]:
            super_stab.discard(coord)
        # Disable gauge2 which anti-commute with gauge and super-stabilizers include gauge2
        basis2 = "XZ"[basis == "X"]
        for coord2, gauge2 in list(self.gauges[basis2].items()):
            if anti_commute(gauge, gauge2):
                self.gauges[basis2].pop(coord2)
                self.super_stabs[basis2] = [s for s in self.super_stabs[basis2] if coord2 not in s]

    def _check(self):
        flag = True
        while flag:
            flag = False

            for coord in self.data_coords:
                if all(coord in self.gauges[basis].keys() for basis in "XZ"):
                    for basis in "XZ":
                        for gauge in self.gauges[basis].values():
                            gauge.discard(coord)
                    flag = True

            for basis, basis2 in ["XZ", "ZX"]:
                used_gauges = set.union(*self.super_stabs[basis]) if self.super_stabs[basis] else set()
                for coord, gauge in list(self.gauges[basis].items()):
                    if not any(anti_commute(gauge, gauge2) for gauge2 in self.gauges[basis2].values()):
                        # Fix a gauge if it is a stabilizer
                        self._fix_gauge(basis, coord)
                        flag = True
                    elif coord not in used_gauges:
                        # Remove unused gauges
                        self.gauges[basis].pop(coord)
                        flag = True

            for basis in "XZ":
                # Update self.gauges
                for coord, gauge in list(self.gauges[basis].items()):
                    if len(gauge) == 1 and gauge != {coord}:
                        q = gauge.pop()
                        self.gauges[basis][q] = {q}
                        for super_stab in self.super_stabs[basis]:
                            if coord in super_stab:
                                super_stab.symmetric_difference_update({q})
                        flag = True
                    if len(gauge) == 0:
                        self.gauges[basis].pop(coord)
                        for super_stab in self.super_stabs[basis]:
                            super_stab.discard(coord)

                # Update self.stabs
                for coord, stab in list(self.stabs[basis].items()):
                    if len(stab) == 1:
                        q = stab.pop()
                        for measurement in chain(self.stabs[basis].values(), self.gauges[basis].values()):
                            measurement.discard(q)
                        self.observable[basis].discard(q)
                        for k in range(2):
                            self.edges[basis][k].discard(q)
                        flag = True
                    if len(stab) == 0:
                        self.stabs[basis].pop(coord)

                # Update self.super_stabs
                self.super_stabs[basis] = [super_stab for super_stab in self.super_stabs[basis] if super_stab]

            ## the following two parts are just for the special case which occurs when defect size is large. It is time-costing and can be removed when defects are not denese.

            # Split super-stabilizer
            if not flag:
                for basis, basis2 in ["XZ", "ZX"]:
                    anti_comm_table = {coord: set() for coord in self.gauges[basis].keys()}  # {basis: {basis2}}
                    for coord, gauge in self.gauges[basis].items():
                        for coord2, gauge2 in self.gauges[basis2].items():
                            if anti_commute(gauge, gauge2):
                                anti_comm_table[coord].add(coord2)

                    for gauge_coords in self.super_stabs[basis]:
                        coord = gauge_coords.pop()
                        super_stab = {coord}
                        anti_comm_gauges = anti_comm_table[coord].copy()
                        while anti_comm_gauges:
                            loop_flag = False
                            for coord, gauges2 in anti_comm_table.items():
                                if coord not in super_stab and coord in gauge_coords and gauges2.intersection(
                                        anti_comm_gauges):
                                    gauge_coords.remove(coord)
                                    super_stab.add(coord)
                                    anti_comm_gauges.symmetric_difference_update(gauges2)
                                    loop_flag = True
                            assert loop_flag

                        if gauge_coords:
                            self.super_stabs[basis].append(super_stab)
                            flag = True
                        else:
                            gauge_coords.update(super_stab)

            # delete seperate part
            if not flag:
                measurements = {basis:{**self.stabs[basis], **self.gauges[basis]} for basis in "XZ"}
                for basis in "XZ":
                    G = nx.Graph()
                    G.add_nodes_from(measurements[basis].keys())

                    G_edges = {q: set() for q in self.data_coords}
                    for coord, measurement in measurements[basis].items():
                        for q in measurement:
                            G_edges[q].add(coord)
                    for q in self.data_coords:
                        if len(G_edges[q]) == 2:
                            G.add_edge(*G_edges[q])

                    for idx, gauge_coords in enumerate(self.super_stabs[basis]):
                        G.add_node(idx)
                        for coord in gauge_coords:
                            G.add_edge(idx, coord)

                    for component in list(nx.connected_components(G)):
                        data_qubits_component = set().union(*(measurements[basis][coord] for coord in component if type(coord) is tuple))

                        super_stab = set()
                        for coord in component:
                            if type(coord) is tuple:
                                super_stab.symmetric_difference_update(measurements[basis][coord])

                        if any(super_stab.issubset(self.edges[basis][k]) for k in range(2)):
                            # delete data_qubits_component
                            for basis2 in "XZ":
                                for measurement in measurements[basis2].values():
                                    measurement.difference_update(data_qubits_component)
                                for k in range(2):
                                    self.edges[basis2][k].difference_update(data_qubits_component)
                                self.observable[basis2].difference_update(data_qubits_component)
                            flag = True
                        # else:
                        #     assert all(self.corners[i][j] in data_qubits_component for i, j in product(range(2), repeat=2))

                self.data_coords = set().union(*(m for basis in "XZ" for m in chain(self.stabs[basis].values(),
                                                                                    self.gauges[basis].values(),
                                                                                    [self.observable[basis]],
                                                                                    )))

                self.qubit_coords = self.data_coords.union(*(set(chain(self.stabs[basis].keys(),
                                                                       self.gauges[basis].keys())) for basis in "XZ"))
        for i, j in product(range(2), repeat=2):
            assert not self.edges["X"][i].intersection(self.edges["Z"][j]).difference({self.corners[i][j]})

if __name__ == "__main__":
    d = 15
    defects = [(20, 20), (3, 13), (28, 12),
            (19, 9), (10, 6), (5, 19), (8, 18), (17, 21), (11, 23), (13, 17), (21, 9),
           (15, 23), (24, 26), (16, 22), (22, 10), (5, 3), (8, 2), (3, 15), (28, 14),
            (17, 23), (2, 4), (0, 16),
           (13, 1),
           (26, 16), (14, 8), (5, 5), (9, 3),
           (3, 17), (28, 16), (23, 29), (9, 21), (15, 9), (6, 6), (1, 19), (26, 18),
           (18, 14), (25, 29), (16, 26), (22, 14), (29, 29), (5, 7), (20, 26), (21, 25),
           (12, 22), (4, 18), (13, 5), (26, 2), (24, 14),
           (13, 23),
           (18, 16), (29, 13),
           (12, 6), (28, 2), (22, 16), (14, 12), (17, 11), (3, 21), (10, 8), (1, 5),
           (26, 4), (8, 20), (15, 13), (7, 9), (10, 26), (2, 22), (29, 15), (16, 30),
           (21, 11), (3, 5), (22, 18), (4, 4), (14, 14), (5, 11), (19, 13), (10, 10),
           (1, 7), (13, 9), (8, 22), (18, 2), (25, 17), (7, 11), (1, 25), (23, 1),
           (29, 17), (20, 14),
           (6, 24), (4, 6),

           (27, 29),
           (29, 25),
           (3, 25),
           (29, 23),
            (29, 21),
            (27, 21),
            (25, 21),
            (13, 21),
            (29, 1),
            (23, 3),
            (25, 5),
            (5, 17),
            (7, 15),
            (9, 15),
               (11, 21)
# (11, 17),
            ]
    Q = LogicalQubit(d, True)

    for q in defects:
        Q.disable(q)

    Q.disable((11, 17))
    Q.update_distance()
    print(Q.distance)
    pass
