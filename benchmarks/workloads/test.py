from pathlib import Path

import stim

from tqecd.construction import annotate_detectors_automatically

HERE = Path(__file__).parent
MAIN_BENCHMARKS_FOLDER = HERE.parent
WORKLOADS_FOLDER = MAIN_BENCHMARKS_FOLDER / "workloads"

BENCHMARK_SITUATIONS: dict[str, list[Path]] = {
    folder.name: list(file for file in folder.iterdir())
    for folder in WORKLOADS_FOLDER.iterdir()
    if folder.is_dir()
}
ALL_BENCHMARK_FILES: list[Path] = sum(BENCHMARK_SITUATIONS.values(), start=[])


for filepath in sorted(ALL_BENCHMARK_FILES):
    print(filepath)
    circuit = stim.Circuit.from_file(filepath)
    annotated_circuit = annotate_detectors_automatically(circuit)
