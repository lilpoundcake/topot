"""Tests for src/topot/utils/validator.py — strict input validation.

Positive controls run the full CLI on every existing test case.
Negative controls craft minimal broken inputs in tmp_path and confirm the
validator rejects them with an actionable error.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

# Allow tests to import the package without an install step.
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / "src"))

from topot.utils.validator import (  # noqa: E402
    ValidationError,
    _validate_dual_topology,
    _validate_ff_dir,
    _validate_gro_header,
    _validate_local_itp_includes,
    _validate_ndx,
    _validate_top_structure,
    validate_inputs_pre_parse,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _run_cli(args, expect_exit=0):
    """Invoke `python -m topot.cli` with the source tree on PYTHONPATH."""
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    proc = subprocess.run(
        [sys.executable, "-m", "topot.cli", *map(str, args)],
        capture_output=True, text=True, env=env,
    )
    assert proc.returncode == expect_exit, (
        f"Expected exit {expect_exit} but got {proc.returncode}.\n"
        f"stdout: {proc.stdout}\nstderr: {proc.stderr}"
    )
    return proc


def _make_fake_gro(path, atom_count_declared, atom_records):
    """Write a minimal .gro file with custom declared count and records.

    Each record is a fixed-width GROMACS line with placeholder atom data.
    """
    lines = ["PMX MODEL\n", f"{atom_count_declared}\n"]
    for i in range(1, atom_records + 1):
        lines.append(
            f"{1:>5d}{'ALA':<5}{'CA':>5}{i:>5d}"
            f"{0.0:>8.3f}{0.0:>8.3f}{0.0:>8.3f}\n"
        )
    lines.append("   1.00000   1.00000   1.00000\n")
    path.write_text("".join(lines))


def _make_minimal_top(path, includes, molecules):
    """Write a minimal .top with the requested includes and [ molecules ]."""
    body = []
    for inc in includes:
        body.append(f'#include "{inc}"\n')
    body.append("\n[ system ]\nTest\n\n[ molecules ]\n")
    for name, count in molecules:
        body.append(f"{name} {count}\n")
    path.write_text("".join(body))


def _make_minimal_itp(path, moltype_name, n_atoms, with_typeB=True):
    """Write a minimal .itp with [ moleculetype ] and [ atoms ]."""
    lines = [
        "[ moleculetype ]\n",
        f"{moltype_name} 3\n",
        "\n[ atoms ]\n",
    ]
    for i in range(1, n_atoms + 1):
        if with_typeB:
            lines.append(
                f"{i} N 1 ALA N {i} 0.0 14.0 N 0.0 14.0\n"
            )
        else:
            lines.append(f"{i} N 1 ALA N {i} 0.0 14.0\n")
    path.write_text("".join(lines))


# ---------------------------------------------------------------------------
# Positive controls — every existing test case must still process end-to-end
# ---------------------------------------------------------------------------

TEST_CASES = [
    # (test_dir, gro_name, with_ndx, lambda_0_expected, lambda_1_expected)
    ("H_TRP33TYR",                          "md_mut.gro", True,  6487, 6484),
    ("I_MET30LYS-I_LEU29ARG",               "md.gro",     False, 3966, 3976),
    ("A_ARG155ASH-A_ASP177ASH-A_LYS180ASP", "md.gro",     False, 6282, 6262),
]


@pytest.mark.parametrize("test_dir,gro_name,with_ndx,l0_count,l1_count", TEST_CASES)
def test_cli_passes_validation_on_real_inputs(tmp_path, test_dir, gro_name,
                                              with_ndx, l0_count, l1_count):
    data_dir = ROOT / "tests" / test_dir
    args = ["-g", data_dir / gro_name, "-p", data_dir / "topol.top",
            "-o", tmp_path / "out"]
    if with_ndx:
        args += ["-n", data_dir / "index.ndx"]

    proc = _run_cli(args)
    assert "Step 1: Validating input files..." in proc.stdout
    assert "Step 6:" in proc.stdout
    # Protein-only output should match documented expectations exactly.
    l0_text = (tmp_path / "out" / "md_lambda_0.gro").read_text().splitlines()
    l1_text = (tmp_path / "out" / "md_lambda_1.gro").read_text().splitlines()
    assert int(l0_text[1]) == l0_count
    assert int(l1_text[1]) == l1_count


# ---------------------------------------------------------------------------
# Negative controls — each check should reject its bad input
# ---------------------------------------------------------------------------

def test_gro_atom_count_mismatch(tmp_path):
    gro = tmp_path / "bad.gro"
    _make_fake_gro(gro, atom_count_declared=10, atom_records=5)
    with pytest.raises(ValidationError, match="atom count mismatch"):
        _validate_gro_header(gro)


def test_gro_non_integer_count(tmp_path):
    gro = tmp_path / "bad.gro"
    gro.write_text("title\nnot_a_number\nfoo\n")
    with pytest.raises(ValidationError, match="not an integer"):
        _validate_gro_header(gro)


def test_top_missing_molecules_section(tmp_path):
    top = tmp_path / "bad.top"
    top.write_text('#include "chain.itp"\n[ system ]\nT\n')
    with pytest.raises(ValidationError, match=r"\[ molecules \] section"):
        _validate_top_structure(top)


def test_top_no_includes_and_no_moltype(tmp_path):
    top = tmp_path / "bad.top"
    top.write_text("[ system ]\nT\n[ molecules ]\nProtein 1\n")
    with pytest.raises(ValidationError, match="no #include"):
        _validate_top_structure(top)


def test_dangling_local_itp_include(tmp_path):
    top = tmp_path / "topol.top"
    _make_minimal_top(top, includes=["missing_chain.itp"],
                      molecules=[("Protein_chain_A", 1)])
    with pytest.raises(ValidationError, match="Missing local #include"):
        _validate_local_itp_includes(top, ["missing_chain.itp"])


def test_ff_includes_are_not_required_locally(tmp_path):
    """Force-field includes (.ff/ paths, tip3p.itp, ions.itp) must not raise."""
    top = tmp_path / "topol.top"
    top.write_text("")
    # All of these should be silently skipped.
    _validate_local_itp_includes(top, [
        "amber99sb-star-ildn-mut.ff/forcefield.itp",
        "./amber99sb-star-ildn-mut.ff/tip3p.itp",
        "./amber99sb-star-ildn-mut.ff/ions.itp",
        "../shared/forcefield.itp",
        "posre_chain_A.itp",
    ])


def test_dual_topology_rejects_non_pmx_file():
    atoms = {1: {'type': 'N', 'typeB': None},
             2: {'type': 'C', 'typeB': None}}
    with pytest.raises(ValidationError, match="Not a dual-topology file"):
        _validate_dual_topology(atoms)


def test_dual_topology_accepts_pmx_file():
    atoms = {1: {'type': 'N', 'typeB': 'N'},
             2: {'type': 'DUM_C', 'typeB': 'C'}}
    _validate_dual_topology(atoms)  # must not raise


def test_ff_dir_does_not_exist(tmp_path):
    with pytest.raises(ValidationError, match="not found"):
        _validate_ff_dir(tmp_path / "does_not_exist")


def test_ff_dir_empty(tmp_path):
    empty = tmp_path / "empty_ff"
    empty.mkdir()
    with pytest.raises(ValidationError, match="no force fields"):
        _validate_ff_dir(empty)


def test_ff_dir_with_ff_subdir(tmp_path):
    parent = tmp_path / "mutff"
    (parent / "fake.ff").mkdir(parents=True)
    _validate_ff_dir(parent)  # must not raise


def test_ndx_without_groups(tmp_path):
    ndx = tmp_path / "bad.ndx"
    ndx.write_text("# just a comment, no groups here\n1 2 3 4 5\n")
    with pytest.raises(ValidationError, match="no group headers"):
        _validate_ndx(ndx)


def test_ndx_with_groups_passes(tmp_path):
    ndx = tmp_path / "good.ndx"
    ndx.write_text("[ System ]\n1 2 3 4 5\n")
    _validate_ndx(ndx)  # must not raise


# ---------------------------------------------------------------------------
# Integration negative cases — full CLI must exit 1 with a useful message
# ---------------------------------------------------------------------------

def test_cli_rejects_truncated_gro(tmp_path):
    """Truncated GRO: declared count > actual atom records."""
    data_dir = ROOT / "tests" / "H_TRP33TYR"
    bad_gro = tmp_path / "bad.gro"

    src_lines = (data_dir / "md_mut.gro").read_text().splitlines(keepends=True)
    # Keep title + count line, but drop most atom records (and the box-vectors).
    truncated = src_lines[:2] + src_lines[2:12]  # 10 records << declared count
    bad_gro.write_text("".join(truncated))

    proc = _run_cli(
        ["-g", bad_gro, "-p", data_dir / "topol.top", "-o", tmp_path / "out"],
        expect_exit=1,
    )
    assert "atom count mismatch" in proc.stdout


def test_cli_rejects_top_without_molecules(tmp_path):
    data_dir = ROOT / "tests" / "H_TRP33TYR"
    bad_top = tmp_path / "bad.top"
    # Copy itp files alongside the bad top so the include check doesn't
    # short-circuit first.
    for itp in data_dir.glob("*.itp"):
        shutil.copy(itp, tmp_path / itp.name)
    bad_top.write_text(
        '#include "pmx_topol_Protein_chain_B.itp"\n[ system ]\nT\n'
    )
    proc = _run_cli(
        ["-g", data_dir / "md_mut.gro", "-p", bad_top, "-o", tmp_path / "out"],
        expect_exit=1,
    )
    assert "[ molecules ]" in proc.stdout


def test_cli_rejects_missing_itp_include(tmp_path):
    data_dir = ROOT / "tests" / "H_TRP33TYR"
    bad_top = tmp_path / "bad.top"
    _make_minimal_top(bad_top,
                      includes=["pmx_topol_Protein_chain_DOES_NOT_EXIST.itp"],
                      molecules=[("Protein_chain_DOES_NOT_EXIST", 1)])
    proc = _run_cli(
        ["-g", data_dir / "md_mut.gro", "-p", bad_top, "-o", tmp_path / "out"],
        expect_exit=1,
    )
    assert "Missing local #include" in proc.stdout


def test_cli_rejects_non_dual_topology(tmp_path):
    """Topology with [ atoms ] but no typeB columns — not a PMX file."""
    chain_itp = tmp_path / "chain.itp"
    _make_minimal_itp(chain_itp, "Protein_chain_A", n_atoms=3, with_typeB=False)
    top = tmp_path / "topol.top"
    _make_minimal_top(top, includes=["chain.itp"],
                      molecules=[("Protein_chain_A", 1)])

    gro = tmp_path / "md.gro"
    _make_fake_gro(gro, atom_count_declared=3, atom_records=3)

    proc = _run_cli(
        ["-g", gro, "-p", top, "-o", tmp_path / "out"],
        expect_exit=1,
    )
    assert "Not a dual-topology file" in proc.stdout


def test_cli_rejects_empty_ff_dir(tmp_path):
    data_dir = ROOT / "tests" / "H_TRP33TYR"
    empty_ff = tmp_path / "empty_ff"
    empty_ff.mkdir()
    proc = _run_cli(
        ["-g", data_dir / "md_mut.gro", "-p", data_dir / "topol.top",
         "-o", tmp_path / "out", "--ff-dir", empty_ff],
        expect_exit=1,
    )
    assert "no force fields" in proc.stdout


def test_cli_rejects_ndx_without_groups(tmp_path):
    data_dir = ROOT / "tests" / "H_TRP33TYR"
    bad_ndx = tmp_path / "bad.ndx"
    bad_ndx.write_text("1 2 3 4 5\n")
    proc = _run_cli(
        ["-g", data_dir / "md_mut.gro", "-p", data_dir / "topol.top",
         "-n", bad_ndx, "-o", tmp_path / "out"],
        expect_exit=1,
    )
    assert "no group headers" in proc.stdout
