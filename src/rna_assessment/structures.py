from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio.PDB import PDBParser


@dataclass(frozen=True)
class ResidueRecord:
    chain: str
    pos: int
    nt: str
    residue: object

    def key(self) -> str:
        return f"{self.chain}:{self.pos}"


class PDBStructure:
    def __init__(self) -> None:
        self._pdb_file: Path | None = None
        self._structure = None
        self._residues: list[ResidueRecord] = []
        self._selected_indices: list[int] = []
        self._residue_index: dict[str, list[int | None]] = {}
        self._interaction_cache: dict[str, list[tuple[str, int, int, str]]] = {}

    @classmethod
    def from_file(cls, pdb_file: str | Path, index_name: str | Path | None = None) -> "PDBStructure":
        structure = cls()
        structure.load(pdb_file, index_name=index_name)
        return structure

    def load(self, pdb_file: str | Path, index_name: str | Path | None = None) -> bool:
        self._pdb_file = Path(pdb_file)
        parser = PDBParser(QUIET=True)
        self._structure = parser.get_structure("structure", str(self._pdb_file))
        self._residues = []
        self._selected_indices = []
        self._residue_index = {}
        self._interaction_cache = {}

        model = self._structure[0]
        for index, chain in enumerate(model.child_list):
            del index
            for residue in chain.child_list:
                record = ResidueRecord(
                    chain=chain.id,
                    pos=residue.id[1],
                    nt=residue.resname.strip(),
                    residue=residue,
                )
                self._residues.append(record)
                self._residue_index[record.key()] = [len(self._residues) - 1, None]

        if index_name is None:
            self._select_all()
        else:
            self._load_index(Path(index_name))
        return True

    @property
    def pdb_file(self) -> str:
        if self._pdb_file is None:
            raise ValueError("Structure has not been loaded.")
        return str(self._pdb_file)

    @property
    def struct(self):
        if self._structure is None:
            raise ValueError("Structure has not been loaded.")
        return self._structure

    @property
    def res_seq(self) -> list[int]:
        return list(self._selected_indices)

    @property
    def res_list(self) -> list[ResidueRecord]:
        return list(self._residues)

    def raw_sequence(self) -> str:
        return "".join(self._residues[index].nt for index in self._selected_indices)

    def res_sequence(self) -> list[object]:
        return [self._residues[index].residue for index in self._selected_indices]

    def rank_of(self, chain: str, pos: int) -> int | None:
        entry = self._residue_index.get(f"{chain}:{pos}")
        return None if entry is None else entry[1]

    def cached_interactions(self, cache_key: str) -> list[tuple[str, int, int, str]] | None:
        return self._interaction_cache.get(cache_key)

    def set_cached_interactions(
        self, cache_key: str, interactions: list[tuple[str, int, int, str]]
    ) -> None:
        self._interaction_cache[cache_key] = interactions

    def _select_all(self) -> None:
        for index, residue in enumerate(self._residues):
            self._selected_indices.append(index)
            self._residue_index[residue.key()][1] = len(self._selected_indices) - 1

    def _load_index(self, index_file: Path) -> None:
        entries: list[tuple[str, int, int]] = []
        for row in index_file.read_text(encoding="utf-8").splitlines():
            stripped = row.strip()
            if not stripped or stripped.startswith("#"):
                continue
            for chunk in stripped.split(","):
                chain, pos, count = chunk.split(":")
                entries.append((chain, int(pos), int(count)))

        self._selected_indices = []
        for chain, pos, count in entries:
            start_index = self._absolute_index(chain, pos)
            if start_index is None:
                raise KeyError(f"Bad index key: '{chain}:{pos}'.")
            for absolute_index in range(start_index, start_index + count):
                if absolute_index >= len(self._residues):
                    raise IndexError(f"Index count overflows structure for '{chain}:{pos}:{count}'.")
                record = self._residues[absolute_index]
                if record.chain != chain:
                    raise IndexError(
                        f"Index '{chain}:{pos}:{count}' crosses into chain '{record.chain}'."
                    )
                self._selected_indices.append(absolute_index)
                self._residue_index[record.key()][1] = len(self._selected_indices) - 1

    def _absolute_index(self, chain: str, pos: int) -> int | None:
        entry = self._residue_index.get(f"{chain}:{pos}")
        return None if entry is None else int(entry[0])
