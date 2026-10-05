from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Sequence, Tuple


@dataclass(frozen=True)
class AggregatedJSONSource:
    """One aggregated PCMDI ENSO JSON dataset and its observation choice."""

    data_dir: Path
    dataset: str
    observation: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "data_dir", Path(self.data_dir))
        if len(self.dataset.split(".")) != 3:
            raise ValueError(
                "dataset must have the form '<mip>.<experiment>.<case_id>', "
                f"got {self.dataset!r}"
            )
        if not self.observation:
            raise ValueError("observation must not be empty")

    @property
    def components(self) -> Tuple[str, str, str]:
        mip, experiment, case_id = self.dataset.split(".")
        return mip, experiment, case_id

    def collection_dir(self, collection: str) -> Path:
        return self.data_dir.joinpath(*self.components, collection)

    def resolve_latest_files(self, collections: Sequence[str]) -> Dict[str, str]:
        resolved = {}
        for collection in collections:
            directory = self.collection_dir(collection)
            files = sorted(
                directory.glob("*.json"),
                key=lambda path: (path.stat().st_mtime, path.name),
            )
            if not files:
                raise FileNotFoundError(f"No JSON files found under {directory}")
            resolved[collection] = str(files[-1])
        return resolved
