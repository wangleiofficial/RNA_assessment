from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

from .exceptions import ToolExecutionError, ToolResolutionError
from .structures import prepare_external_structure_input
from .tools import default_tool_registry


@dataclass(frozen=True)
class MolProbityResult:
    clashscore: float | None
    molprobity_score: float | None
    bond_outliers: float | None
    angle_outliers: float | None
    pucker_outliers: float | None
    suite_outliers: float | None
    binary_path: str


class MolProbityRunner:
    _FLOAT_RE = r"([0-9]+(?:\.[0-9]+)?)"
    _PATTERNS = {
        "clashscore": (
            re.compile(rf"(?im)^\s*(?:all-atom\s+)?clashscore\b[^0-9-]*{_FLOAT_RE}\s*$"),
            re.compile(rf"(?im)\bclashscore\b\s*[:=]\s*{_FLOAT_RE}"),
        ),
        "molprobity_score": (
            re.compile(rf"(?im)^\s*molprobity\s+score\b[^0-9-]*{_FLOAT_RE}\s*$"),
            re.compile(rf"(?im)\bmolprobity\s+score\b\s*[:=]\s*{_FLOAT_RE}"),
        ),
        "bond_outliers": (
            re.compile(rf"(?im)\b(?:bad\s+bonds|bond\s+outliers)\b[^0-9-]*{_FLOAT_RE}"),
        ),
        "angle_outliers": (
            re.compile(rf"(?im)\b(?:bad\s+angles|angle\s+outliers)\b[^0-9-]*{_FLOAT_RE}"),
        ),
        "pucker_outliers": (
            re.compile(rf"(?im)\b(?:pucker\s+outliers|bad\s+puckers)\b[^0-9-]*{_FLOAT_RE}"),
        ),
        "suite_outliers": (
            re.compile(rf"(?im)\b(?:suite\s+outliers|bad\s+suites)\b[^0-9-]*{_FLOAT_RE}"),
        ),
    }

    def __init__(self, binary_path: str | Path | None = None) -> None:
        self.binary_path = Path(binary_path) if binary_path else None

    def resolve_binary(self) -> Path:
        registry = default_tool_registry()
        status = registry.status("molprobity", override=self.binary_path)
        if status.binary_path is not None:
            return Path(status.binary_path)
        raise ToolResolutionError(
            "MolProbity is not available. Set 'RNA_KIT_MOLPROBITY' or pass '--molprobity'."
        )

    def validate(self, structure_file: str | Path) -> MolProbityResult:
        binary = self.resolve_binary()
        source_path = Path(structure_file)

        try:
            if source_path.suffix.lower() in {".cif", ".mmcif"}:
                with TemporaryDirectory(prefix="rna-kit-molprobity-") as temp_dir:
                    prepared_input = prepare_external_structure_input(
                        source_path,
                        Path(temp_dir) / f"{source_path.stem}.pdb",
                    )
                    completed = subprocess.run(
                        [str(binary), str(prepared_input)],
                        check=True,
                        capture_output=True,
                        text=True,
                    )
            else:
                completed = subprocess.run(
                    [str(binary), str(source_path)],
                    check=True,
                    capture_output=True,
                    text=True,
                )
        except OSError as exc:
            raise ToolExecutionError(f"Failed to execute MolProbity at '{binary}'.") from exc
        except subprocess.CalledProcessError as exc:
            message = exc.stderr.strip() or exc.stdout.strip() or "unknown error"
            raise ToolExecutionError(f"MolProbity failed for '{source_path}': {message}") from exc

        return self.parse(completed.stdout, completed.stderr, binary_path=binary)

    def parse(
        self,
        stdout: str,
        stderr: str = "",
        *,
        binary_path: str | Path = "molprobity",
    ) -> MolProbityResult:
        combined = "\n".join(part for part in (stdout.strip(), stderr.strip()) if part).strip()
        if not combined:
            raise ToolExecutionError("MolProbity produced no output.")

        stripped = combined.strip()
        if re.fullmatch(self._FLOAT_RE, stripped):
            return MolProbityResult(
                clashscore=float(stripped),
                molprobity_score=None,
                bond_outliers=None,
                angle_outliers=None,
                pucker_outliers=None,
                suite_outliers=None,
                binary_path=str(binary_path),
            )

        values = {
            key: _first_match(patterns, combined)
            for key, patterns in self._PATTERNS.items()
        }
        if all(value is None for value in values.values()):
            raise ToolExecutionError("MolProbity output could not be parsed into summary metrics.")

        return MolProbityResult(
            clashscore=values["clashscore"],
            molprobity_score=values["molprobity_score"],
            bond_outliers=values["bond_outliers"],
            angle_outliers=values["angle_outliers"],
            pucker_outliers=values["pucker_outliers"],
            suite_outliers=values["suite_outliers"],
            binary_path=str(binary_path),
        )


def calculate_molprobity(
    structure_file: str | Path,
    runner: MolProbityRunner | None = None,
) -> MolProbityResult:
    runner = runner or MolProbityRunner()
    return runner.validate(structure_file)


def _first_match(patterns: tuple[re.Pattern[str], ...], text: str) -> float | None:
    for pattern in patterns:
        match = pattern.search(text)
        if match is not None:
            return float(match.group(1))
    return None
