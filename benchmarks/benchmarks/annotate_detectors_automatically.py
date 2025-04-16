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


class BaseAutomaticDetectorAnnotation:
    params = ALL_BENCHMARK_FILES

    def setup(self, path: Path):
        self.circuit = stim.Circuit.from_file(path)

    def teardown(self, _: Path):
        del self.circuit

    def _annotate(self) -> None:
        annotate_detectors_automatically(self.circuit)


class TimeAutomaticDetectorAnnotation(BaseAutomaticDetectorAnnotation):
    def time_annotate_detectors_automatically(self, _: Path):
        self._annotate()


class MemAutomaticDetectorAnnotation(BaseAutomaticDetectorAnnotation):
    def mem_annotate_detectors_automatically(self, _: Path):
        self._annotate()


class PeakMemAutomaticDetectorAnnotation(BaseAutomaticDetectorAnnotation):
    def peakmem_annotate_detectors_automatically(self, _: Path):
        self._annotate()
