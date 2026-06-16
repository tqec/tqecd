from tqecd.boundary import BoundaryStabilizer
from tqecd.flow import FragmentFlows, FragmentLoopFlows
from tqecd.match import _shallow_copy_flows
from tqecd.pauli import PauliString


def _flows() -> FragmentFlows:
    X0 = PauliString({0: "X"})
    Z1 = PauliString({1: "Z"})
    creation = [BoundaryStabilizer(X0, [], [], frozenset([0]), True)]
    destruction = [BoundaryStabilizer(Z1, [], [], frozenset([1]), False)]
    return FragmentFlows(creation, destruction, total_number_of_measurements=2)


def test_shallow_copy_flows_duplicates_lists_but_shares_stabilizers() -> None:
    flows = _flows()
    copied = _shallow_copy_flows(flows)
    assert isinstance(copied, FragmentFlows)

    # The lists are independent: mutating the copy leaves the original intact.
    assert copied.creation is not flows.creation
    assert copied.destruction is not flows.destruction
    copied.creation.pop()
    assert len(flows.creation) == 1
    assert copied.total_number_of_measurements == flows.total_number_of_measurements

    # The boundary stabilizers themselves are shared (never mutated in place).
    other = _shallow_copy_flows(flows)
    assert isinstance(other, FragmentFlows)
    assert other.creation[0] is flows.creation[0]


def test_shallow_copy_flows_recurses_into_loops() -> None:
    inner = _flows()
    loop = FragmentLoopFlows([inner], repeat=3)
    copied = _shallow_copy_flows(loop)
    assert isinstance(copied, FragmentLoopFlows)
    assert copied.repeat == 3
    assert copied.fragment_flows is not loop.fragment_flows
    assert isinstance(copied.fragment_flows[0], FragmentFlows)
    assert copied.fragment_flows[0].creation is not inner.creation
    assert copied.fragment_flows[0].creation[0] is inner.creation[0]
