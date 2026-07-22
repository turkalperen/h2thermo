"""Generate a hydrogen-air property table over the full operating envelope."""

import time

from h2thermo import GridSpecification, ThermoTable, report_progress


def main() -> None:
    """Generate, save and summarise a production-scale table."""
    grid = GridSpecification.linear(
        temperature_range=(200.0, 3000.0),
        pressure_range=(1.0e5, 60.0e5),
        equivalence_ratio_range=(0.2, 1.0),
        shape=(50, 20, 20),
    )
    print(f"Grid shape {grid.shape}, {grid.size} nodes")

    start = time.perf_counter()
    table = ThermoTable.generate(grid, progress=report_progress)
    elapsed = time.perf_counter() - start

    path = table.save("data/generated/h2_air_table.npz")
    size_mb = path.stat().st_size / 1.0e6

    print(f"Generated in {elapsed:.1f} s")
    print(f"Non-convergent nodes: {table.failed_node_count}")
    print(f"Saved {path} ({size_mb:.2f} MB)")


if __name__ == "__main__":
    main()