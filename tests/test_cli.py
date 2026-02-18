"""Integration tests for topot CLI with H_TRP33TYR mutation"""

import tempfile
from pathlib import Path
import sys
import subprocess

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


def test_cli_with_h_trp33tyr_data():
    """Test CLI with H_TRP33TYR test data (W2Y mutation at position 337 in chain H)"""
    test_data_dir = Path(__file__).parent / "H_TRP33TYR"
    mutff_dir = Path(__file__).parent.parent / "mutff"

    # Check if test data exists
    assert test_data_dir.exists(), f"Test data directory not found: {test_data_dir}"

    # Accept both md.gro and md_mut.gro
    gro_file = test_data_dir / "md_mut.gro"
    if not gro_file.exists():
        gro_file = test_data_dir / "md.gro"

    assert gro_file.exists(), f"GRO file not found: {gro_file}"
    assert (test_data_dir / "newtop.top").exists(), "Topology file not found"
    assert (test_data_dir / "index.ndx").exists(), "Index file not found"

    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir)

        # Build command
        cmd = [
            "topot",
            "-g", str(gro_file),
            "-p", str(test_data_dir / "newtop.top"),
            "-n", str(test_data_dir / "index.ndx"),
            "-o", str(output_dir),
        ]

        # Add FF directory if it exists
        if mutff_dir.exists():
            cmd.extend(["--ff-dir", str(mutff_dir)])

        # Run CLI via subprocess
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"CLI returned non-zero exit code: {result.stderr}"

        # Check output files exist
        assert (output_dir / "md_lambda_0.gro").exists(), "md_lambda_0.gro not found"
        assert (output_dir / "md_lambda_1.gro").exists(), "md_lambda_1.gro not found"
        assert (output_dir / "md_lambda_0_WI.gro").exists(), "md_lambda_0_WI.gro not found"
        assert (output_dir / "md_lambda_1_WI.gro").exists(), "md_lambda_1_WI.gro not found"
        assert (output_dir / "md_lambda_0.pdb").exists(), "md_lambda_0.pdb not found"
        assert (output_dir / "md_lambda_1.pdb").exists(), "md_lambda_1.pdb not found"
        assert (output_dir / "index.ndx").exists(), "index.ndx not found"

        # Check atom counts (check _WI.gro which includes water/ions)
        with open(output_dir / "md_lambda_0_WI.gro") as f:
            lines = f.readlines()
            lambda_0_count = int(lines[1].strip())

        with open(output_dir / "md_lambda_1_WI.gro") as f:
            lines = f.readlines()
            lambda_1_count = int(lines[1].strip())

        # For H_TRP33TYR (W2Y mutation):
        # Total atoms: 83,949
        # Lambda 0: ~83,938 (11 DUM_ atoms removed)
        # Lambda 1: ~83,935 (14 DUM_ atoms removed)
        assert lambda_0_count >= 83900, f"Lambda 0 atom count too low: {lambda_0_count}"
        assert lambda_1_count >= 83900, f"Lambda 1 atom count too low: {lambda_1_count}"
        assert lambda_0_count > lambda_1_count, f"Lambda 0 ({lambda_0_count}) should have more atoms than Lambda 1 ({lambda_1_count})"

        # Check protein counts in PDB
        with open(output_dir / "md_lambda_0.pdb") as f:
            l0_pdb_count = sum(1 for line in f if line.startswith("ATOM"))

        with open(output_dir / "md_lambda_1.pdb") as f:
            l1_pdb_count = sum(1 for line in f if line.startswith("ATOM"))

        # Expected: ~6487 protein atoms at lambda 0, ~6484 at lambda 1
        assert l0_pdb_count > 6400, f"Lambda 0 protein atoms too low: {l0_pdb_count}"
        assert l1_pdb_count > 6400, f"Lambda 1 protein atoms too low: {l1_pdb_count}"

        # Check index file has lambda groups
        with open(output_dir / "index.ndx") as f:
            ndx_content = f.read()
            assert "lambda_0" in ndx_content, "lambda_0 group not found in index file"
            assert "lambda_1" in ndx_content, "lambda_1 group not found in index file"

        print(f"✓ GRO file: {gro_file.name}")
        print(f"✓ Lambda 0: {lambda_0_count} atoms (total), {l0_pdb_count} atoms (protein)")
        print(f"✓ Lambda 1: {lambda_1_count} atoms (total), {l1_pdb_count} atoms (protein)")
        print(f"✓ Output directory: {output_dir}")
        print("✓ All tests passed!")


if __name__ == "__main__":
    test_cli_with_h_trp33tyr_data()
