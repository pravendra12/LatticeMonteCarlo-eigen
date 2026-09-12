#!/usr/bin/env python3
"""Build and run deterministic tests of the production KMC tracer implementation.

Usage: python3 script/test_tracer_msd.py [configured-cmake-build-directory]
Requires the project's compiler/dependencies and a CMake compile_commands.json.
A predictor stand-in avoids loading a fitted model: this harness supplies events.
All compilation and simulation output goes into a temporary directory.
"""
import csv
import gzip
import json
import math
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
build = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT / "build"
commands = json.loads((build / "compile_commands.json").read_text())
entry = next(c for c in commands if c["file"].endswith("KineticMcAbstract.cpp"))
args = shlex.split(entry["command"])
includes = []
i = 0
while i < len(args):
    if args[i] in ("-I", "-isystem"):
        includes.extend(args[i:i + 2])
        i += 2
    else:
        if args[i].startswith("-I"):
            includes.append(args[i])
        i += 1

with tempfile.TemporaryDirectory(prefix="lmc-tracer-msd-") as temp:
    work = Path(temp)
    sources = ["lmc/testing/tracer_msd.cpp", "lmc/mc/src/KineticMcAbstract.cpp",
               "lmc/mc/src/McAbstract.cpp", "lmc/mc/src/JumpEvent.cpp",
               "lmc/mc/src/ThermodynamicAveraging.cpp", "lmc/cfg/src/Config.cpp",
               "lmc/pred/src/TimeTemperatureInterpolator.cpp"]
    objects = []
    for source in sources:
        obj = work / (Path(source).stem + ".o")
        objects.append(str(obj))
        subprocess.run([args[0], "-std=c++20", "-O0", "-fopenmp",
                        "-ffunction-sections", "-fdata-sections",
                        "-I" + str(ROOT / "lmc/testing/tracer_stubs"), *includes,
                        "-c", str(ROOT / source), "-o", str(obj)], check=True)
    exe = work / "test"
    subprocess.run([args[0], "-fopenmp", *objects, "-Wl,--gc-sections",
                    "-lboost_iostreams", "-lboost_filesystem", "-lboost_system",
                    "-lmpi_cxx", "-lmpi", "-o", str(exe)], check=True)

    def run(mode, directory=None, mpi=False):
        """Run one MPI lifetime, optionally appending to an existing run."""
        directory = directory or work / mode
        directory.mkdir(exist_ok=True)
        command = [str(exe), mode]
        if mpi:
            command = ["mpirun", "-np", "2", *command]
        subprocess.run(command, cwd=directory, check=True)
        return directory

    def check_log(directory, steps):
        """Compare every scalar MSD to the known backtracking trajectory."""
        with (directory / "kmc_log.txt").open() as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        assert [int(r["steps"]) for r in rows] == steps
        assert list(rows[0])[-4:] == ["dR_Fe", "dR_Ni", "MSD_Fe", "MSD_Ni"]
        for row in rows:
            step = int(row["steps"])
            assert abs(float(row["MSD_Fe"]) - 2 * (step % 2)) < 1e-12
            assert float(row["MSD_Ni"]) == 0
            assert float(row["time"]) == step * 0.25
        for path in directory.glob("*.displacements.gz"):
            with gzip.open(path, "rt") as stream:
                atoms = list(csv.DictReader(
                    (line for line in stream if not line.startswith("#")), delimiter=" "))
            expected = sum(sum(float(a[c]) ** 2 for c in ("dx", "dy", "dz"))
                           for a in atoms if a["element"] == "Fe") / 2
            assert abs(expected - 2 * (int(path.name.split(".")[0]) % 2)) < 1e-12

    run("tracking")
    baseline = run("baseline")
    check_log(baseline, list(range(8)))
    assert {p.name for p in baseline.glob("*.displacements.gz")} == {
        "0.displacements.gz", "3.displacements.gz", "6.displacements.gz", "7.displacements.gz"}
    checkpoint = run("checkpoint")
    run("restart", checkpoint)
    check_log(checkpoint, [0, 1, 2, 3, 3, 4, 5, 6, 7])
    continuous = [float(x) for x in (baseline / "state.txt").read_text().split()]
    resumed = [float(x) for x in (checkpoint / "state.txt").read_text().split()]
    assert len(continuous) == len(resumed)
    assert all(math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12)
               for a, b in zip(continuous, resumed))
    fresh = work / "new_origin"
    fresh.mkdir()
    shutil.copy(checkpoint / "3.cfg.gz", fresh)
    run("new_origin", fresh)
    assert (fresh / "state.txt").read_text().splitlines()[0] == "3 0.75 3 0.75"
    bad = work / "bad_header"
    bad.mkdir()
    for name in ("3.cfg.gz", "3.displacements.gz"):
        shutil.copy(checkpoint / name, bad)
    (bad / "kmc_log.txt").write_text("steps\ttime\n")
    run("bad_header", bad)
    assert (bad / "kmc_log.txt").read_text() == "steps\ttime\n"
    adaptive = run("adaptive")
    check_log(adaptive, list(range(8)))
    mpi = run("baseline", work / "mpi", mpi=True)
    check_log(mpi, list(range(8)))
    print("PASS: tracer hops, backtracking, multispecies brute-force checks, restart, "
          "new origin, log columns, checkpoint cadence, legacy header, adaptive and two MPI ranks")
