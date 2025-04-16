from pathlib import Path

import stim

HERE = Path(__file__).parent
CODE_TASKS = [
    "repetition_code:memory",
    "surface_code:rotated_memory_x",
    "surface_code:rotated_memory_z",
    "surface_code:unrotated_memory_x",
    "surface_code:unrotated_memory_z",
    "color_code:memory_xyz",
]
DISTANCES = [3, 5, 7, 9]


def remove_annotations(circuit: stim.Circuit) -> stim.Circuit:
    ret = stim.Circuit()
    for inst in circuit:
        if isinstance(inst, stim.CircuitRepeatBlock):
            ret.append(
                stim.CircuitRepeatBlock(
                    inst.repeat_count, remove_annotations(inst.body_copy())
                )
            )
        elif inst.name not in ["DETECTOR", "OBSERVABLE_INCLUDE", "SHIFT_COORDS"]:
            ret.append(inst)
    return ret


def main() -> None:
    for task in CODE_TASKS:
        folder = HERE / task.replace(":", "_")
        if not folder.exists():
            folder.mkdir(parents=True)

        for distance in DISTANCES:
            circuit = stim.Circuit.generated(task, distance=distance, rounds=distance)
            circuit = remove_annotations(circuit)
            with open(folder / f"{distance}.stim", "w") as f:
                f.write(str(circuit))
    print(
        f"Circuit generated in {HERE}. Hand-modifications are still needed "
        "(replacement of combined instructions by split M and R instructions)."
    )


if __name__ == "__main__":
    main()
