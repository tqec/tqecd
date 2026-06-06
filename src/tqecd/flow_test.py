from copy import copy

import stim

from tqecd.flow import (
    FragmentFlows,
    FragmentLoopFlows,
    _build_flows_from_fragment,
)
from tqecd.fragment import Fragment
from tqecd.measurement import RelativeMeasurementLocation


def test_repeated_measurements_keep_first_offset_and_order() -> None:
    fragment = Fragment(stim.Circuit("R 0\nTICK\nM 0\nTICK\nM 0"))

    flows = _build_flows_from_fragment(fragment)

    first_measurement = RelativeMeasurementLocation(-2, 0)
    assert flows.creation[0].measurements == [
        first_measurement,
        first_measurement,
    ]
    assert [flow.measurements for flow in flows.destruction] == [
        [first_measurement],
        [first_measurement],
    ]


def test_copying_flows_creates_independent_flow_lists() -> None:
    fragment = Fragment(stim.Circuit("R 0\nTICK\nM 0"))
    flows = _build_flows_from_fragment(fragment)

    copied_flows = copy(flows)

    assert isinstance(copied_flows, FragmentFlows)
    assert copied_flows is not flows
    assert copied_flows.creation is not flows.creation
    assert copied_flows.destruction is not flows.destruction
    assert copied_flows.creation[0] is flows.creation[0]

    loop_flows = FragmentLoopFlows([flows], repeat=2)
    copied_loop_flows = copy(loop_flows)

    assert copied_loop_flows is not loop_flows
    assert copied_loop_flows.fragment_flows is not loop_flows.fragment_flows
    assert copied_loop_flows.fragment_flows[0] is not flows
    copied_loop_flows.remove_creation(0)
    assert len(flows.creation) == 1
