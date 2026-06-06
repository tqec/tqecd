import stim

from tqecd.flow import _build_flows_from_fragment
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
